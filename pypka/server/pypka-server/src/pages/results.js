import React from "react"

import Layout from "../components/layout"
import HeaderBar from "../components/headerbar"

import LinePlot from "../components/lineplot"

import comp_cluster from "../images/comp_cluster.jpg"

import styles from "./run-pypka.module.css"
import "../styles/thesaas.scss"

import graph_design from "../images/graph_design.png"

import axios from "axios"

import GlobalState from "../context/ThemeContext"

import { ProgressCircle, ProgressLinear } from "../components/progressBar"

function downloadString(text, fileType, fileName) {
  var blob = new Blob([text], { type: fileType });

  var a = document.createElement('a');
  a.download = fileName;
  a.href = URL.createObjectURL(blob);
  a.dataset.downloadurl = [fileType, a.download, a.href].join(':');
  a.style.display = "none";
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  setTimeout(function() { URL.revokeObjectURL(a.href); }, 1500);
}

async function submit_pypka_calculation (object_state) {
  var config = { headers: {  
    'Content-Type': 'application/json',
    'Access-Control-Allow-Origin': '*'}
  }
  const send_json = {
    subID:  object_state.subID,
    pdb:    object_state.pdb,
    inputNamingScheme: object_state.inputNamingScheme,
    outputpKs:          object_state.outputpKs,
    outputfile:         object_state.outputfile,
    outputNamingScheme: object_state.outputNamingScheme,
    outputFilepH:       object_state.outputFilepH,
    pHmin:  object_state.pHmin,
    pHmax:  object_state.pHmax,
    pHstep: object_state.pHstep,
    epsin:  object_state.epsin,
    epsout: object_state.epsout,
    ionic:  object_state.ionic
  }
  console.log(send_json)
  try {
    const response = await axios.post('http://127.0.0.1:5000/submitSim', send_json, config)
    const data = response.data
    console.log(data)
    return data
  } catch (error) {
    console.log(error)
  }
}

async function run_pypka(object, global) {

  const response = await submit_pypka_calculation(object.state)
  
  const titration_curve = response.titration  
  
  object.setState({
    tit_x: titration_curve[0],
    tit_y: titration_curve[1],
    pKas: response.pKas,
    params: response.parameters,
    pdb_out: response.pdb_out,
    titdatarevision: 2
  })

  global.saveSubmission(object.state.subID, response)
  
}

const PKTable = (props) => {
  if (props.pKas) {
    var pKlist = props.pKas.map((sites) =>
    <tr key={sites[2] + sites[1]}>
      <th scope="row">{sites[0]}</th>
      <td>{sites[1]}</td>
      <td>{sites[2]}</td>
      <td>{sites[3]}</td>
    </tr>    
    );
    var waiting = ''
  } else {
    var pKlist = <tr></tr>
    var waiting = <ProgressCircle />
  }
  
  return (
    <div>
      <table className="table table-hover" style={{textAlign: "center"}}>
        <thead>
          <tr>
            <th>Chain</th>
            <th>Residue Name</th>
            <th>Residue Number</th>
            <th>p<em>K</em><sub>a</sub></th>
          </tr>
        </thead>
        <tbody>
          {pKlist}
        </tbody>
      </table>
      {waiting}
    </div>
  )
}


class Results extends React.Component {
    state = this.props.location.state.send_json

    componentDidMount() {
      const global = new GlobalState(this.state.subID)

      console.log(global)

      if (global.state.pKas.length !== 0 && global.state.tit_x.length !== 0 && global.state.tit_y.length !== 0) {
        this.setState({
          pKas: global.state.pKas,
          tit_x: global.state.tit_x,
          tit_y: global.state.tit_y,
          params: global.state.params,
          pdb_out: global.state.pdb_out,
          titdatarevision: 1
        })
      } else {
        run_pypka(this, global)
      }

      var qs,js,q,s,d=document, gi=d.getElementById, ce=d.createElement, gt=d.getElementsByTagName, id="typef_orm_share", b="https://embed.typeform.com/"; 
      if(!gi.call(d,id)){ js=ce.call(d,"script"); js.id=id; js.src=b+"embed.js"; q=gt.call(d,"script")[0]; q.parentNode.insertBefore(js,q) } 
    }

    DownloadTitration = () => {
      var tit_data = 'pH;pKa\n'
      for (const i in this.state.tit_y) {
        const x = this.state.tit_x[i]
        const y = this.state.tit_y[i]
        tit_data = tit_data.concat(`${x};${y}\n`)
      }
      console.log(tit_data)
      downloadString(tit_data, "text/csv", "titration.csv")
    }

    DownloadpKTable = () => {
      var pK_table = 'Chain;Type;Number;pKa\n'
      for (const elem of this.state.pKas) {
        pK_table = pK_table.concat(`${elem}\n`.split(',').join(';'))
      }
      console.log(pK_table)
      downloadString(pK_table, "text/csv", "pKs.csv")
    }

    DownloadParams = () => {
      const params_cleaned = this.state.params.split(',').join('').split(':').join('=').split('{').join('').split('}').join('').split(" '").join("").split("'").join("")
      downloadString(params_cleaned, "text/txt", "params.txt")
    }

    DownloadPDB = () => {      
      downloadString(this.state.pdb_out, "text/txt", `pdb_${this.state.outputFilepH}.out`)
    }
    
    render() {

      console.log(this.state.titdatarevision)
      

      return (                    
  <Layout>
    <HeaderBar image={comp_cluster} 
               title={"Results"} 
               subtitle={""} />
    {this.state.pKas ? '' : <ProgressLinear run_time={this.state.time_estimate} /> }
    

    <section id="basic" className="section bb-1">
        <div className="container">
          <div className="row gap-y">

            <div className="col-md-6">
              <div className={"card shadow-4 " + styles.card}>
                <div className={"card-body " + styles.cardBody}>
                    <h2>{this.state.protein_name}</h2>
                    <div className="row gap-y">
                        <div className="col-md-6">
                            <p>Number of Titrable Sites: {this.state.nsites}</p>
                            <p>Number of Chains: {this.state.nchains}</p>
                            <p>Number of Amino Acids</p>
                            {this.state.pdb_out ? <button type="button" className="btn btn-outline-primary" onClick={this.DownloadPDB}>PDB at pH {this.state.outputFilepH} </button> : ''}
                        </div>
                        <div className="col-md-6">
                            <p>pH range: {this.state.pHmin}-{this.state.pHmax}</p>
                            <p>Protein Dielectric: {this.state.epsin}</p>
                            <p>Solvent Dielectric: {this.state.epsout}</p>
                            {this.state.params ? <button type="button" className="btn btn-outline-primary" onClick={this.DownloadParams}>All Parameters</button>: ''}
                        </div>
                    </div>                    

                </div>
                
                
              </div>
              <br />
              <div className={"card shadow-4 " + styles.card}>
              {this.state.tit_x ? <button type="button" className="btn btn-outline-primary" onClick={this.DownloadTitration} >Download CSV</button> : '' }
                <div className={"card-body " + styles.cardBody}>
                    {this.state.tit_x ? <LinePlot x={this.state.tit_x} y={this.state.tit_y} width={"100px"} revision={this.state.titdatarevision}/> : <div><h5>Titration Curve</h5><ProgressCircle /></div> }
                </div>
              </div>
            </div>            

            <div className="col-md-6">
              <div className={"card shadow-4 " + styles.card}>
                {this.state.tit_x ? <button type="button" className="btn btn-outline-primary" onClick={this.DownloadpKTable} >Download CSV</button> : '' }
                <div className={"card-body " + styles.cardBody}>
                    <PKTable pKas={this.state.pKas} />
                </div>
              </div>
            </div>


          </div>
        </div>
      </section>

      <section className="section">
        <div className="container">
          <div className="row gap-y align-items-center">
            <div className="col-md-6 text-center">
              <img src={graph_design} alt="..." />
            </div>

            <div className="col-md-6 text-center text-md-left" >
              <h2>Help us design better analysis</h2>
              <p className="lead mb-6">Let us know what other plots and data you would like to see</p>
              <p>
                <a className="btn btn-lg btn-round btn-info typeform-share button" href="https://hdcalgarve.typeform.com/to/jGgOaU" data-mode="popup" target="_blank">
                  Reach Out
                </a>
              </p>
            </div>
          </div>
        </div>
      </section>
    
  </Layout>
)}}

export default Results