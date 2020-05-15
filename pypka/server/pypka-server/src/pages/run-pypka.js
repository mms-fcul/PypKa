import { navigate } from "gatsby"
import React from "react"

import Layout from "../components/layout"
import HeaderBar from "../components/headerbar"
import SEO from "../components/seo"

import styles from "./run-pypka.module.css"
import axios from "axios"

import { toast } from 'react-toastify';
import 'react-toastify/dist/ReactToastify.css';

import comp_cluster from "../images/comp_cluster.jpg"

function getTimeEstimate(nsites, pHmax, pHmin, pHstep) {
  return nsites * 1.2 + ((pHmax - pHmin) / pHstep) * 0.5
}

var config = { headers: {  
    'Content-Type': 'application/json',
    'Access-Control-Allow-Origin': '*',
    'Token': Math.random()}
}

async function downloadPDB(pdbid) {
  const response = await axios.get('https://files.rcsb.org/view/' + pdbid + '.pdb')
  return response.data
}

async function getSubID(pdbid) {
  const response = await axios.post('http://127.0.0.1:5000/getSubID', {}, config)
  return response.data.subID
}

async function getNumberOfTitrableSites(pdbfile) {
  try {
    const response = await axios.post('http://127.0.0.1:5000/getTitrableSitesNumber', {PDB: pdbfile}, config)
    console.log(response)
    return response.data
  } catch (error) {
    console.log(error)
    return error
  }
}

async function checkPdbId(object, pdbid) {
  const response = await axios.get('https://www.rcsb.org/pdb/rest/describeMol?structureId=' + pdbid)  

  const parser = new DOMParser();
  const xmlDoc = parser.parseFromString(response.data, "text/xml");
  console.log(xmlDoc)  
  const error = xmlDoc.getElementsByTagName('structureId')[0].getAttribute('error')
  let pdbid_final = pdbid
  let pdb_file = null
  var nchains = 0
  var nsites = 0
  var protein_name = ''
  console.log(error)
  if (error) {
    pdbid_final = null
  } else {
    protein_name = xmlDoc.getElementsByTagName('macroMolecule')[0].getAttribute('name') + ` (${pdbid})`
    pdb_file = await downloadPDB(pdbid)
    var data = await getNumberOfTitrableSites(pdb_file)
                                  
    if ('nchains' in data) {
      nchains = data.nchains
      nsites = data.nsites
    }
    console.log(data, 'nchains' in data)
    console.log(protein_name, nchains, nsites)
  }

  object.setState({
    pdbcode: pdbid_final,
    pdbcode_error: error,
    pdbfile: pdb_file,
    nchains: nchains,
    nsites: nsites,
    protein_name: protein_name
  });
}

async function submit_pypka_calculation (object) {
  object.setTimeEstimate()
  const subID = await getSubID()
  const send_json = {
    subID: subID,
    pdb:               object.state.pdbfile,
    inputNamingScheme: object.state.inputFileNamingScheme,
    outputpKs:         object.state.outputpKValues,
    outputfile:         object.state.outputPDBFile,
    outputNamingScheme: object.state.outputFileNamingScheme,
    outputFilepH:       object.state.outputFilepH,      
    pHmin:   object.state.pHmin,
    pHmax:   object.state.pHmax,
    pHstep:  object.state.pHstep,
    epsin:   object.state.proteinDielectric,
    epsout:  object.state.solventDielectric,
    ionic:   object.state.ionicStrength,
    nchains: object.state.nchains,
    nsites:  object.state.nsites,
    protein_name: object.state.protein_name,
    time_estimate: object.state.time_estimate
  }
  navigate(
    "/results/",
    {
      state: { send_json },
    }
  )
}


const InputFloat = (props) => (
  <div className={styles.inputGroup}> 
    <span className="input-group-addon" style={{width: "50%"}}>{props.name}</span>
    <InputFloatSimple extraStyle={props.extraStyle}
                      defaultValue={props.defaultValue}
                      condition1={props.condition1}
                      condition2={props.condition2}
                      disableOn={props.disableOn}
                      saveState={props.saveState}
                      width={"50%"} />
  </div>
)

const InputFloatSimple = (props) => (
    <input  className={ (props.condition1 && props.condition2) ? styles.borderPrimary + ' ' + props.extraStyle : props.extraStyle }
            style={{width: props.width, textAlign: "center", padding: "5px 12px"}}
            type="text" defaultValue={props.defaultValue} 
            aria-describedby="basic-addon1"
            disabled={!props.disableOn} 
            onChange={(e) => {
              const displayed_value = e.target.value
              const pH_value = parseFloat(displayed_value)
            
              if (displayed_value !== pH_value) {
                if (isNaN(pH_value)) {
                  console.log(pH_value)
                  e.target.value = ''
                } else {
                  console.log(pH_value)
                  e.target.value = pH_value
                }
              }
              props.saveState(pH_value)
            }} />
)


class RunPage extends React.Component {
    state = {
      inputModeRadio: 'pdb_code',
      pdbcode: null,
      pdbcode_error: null,
      inputFile: null,
      inputFileNamingScheme: 'GROMOS',

      outputPDBFile: null,
      outputFileNamingScheme: null,
      outputFilepH: 7,
      outputpKValues: true,

      pHmin: 0,
      pHmax: 12,
      pHstep: 0.25,

      solventDielectric: 80,
      proteinDielectric: 12,

      ionicStrength: 0.1,

      numberOfTitratableSites: 0,
      pdbfile: '',

      nsites: 0,
      nchains: 0,

      time_estimate: 0

    }

    setTimeEstimate = () => {
      var estimate = getTimeEstimate(this.state.nsites, this.state.pHmax, this.state.pHmin, this.state.pHstep)
      this.setState({
        time_estimate: estimate
      })
    }
    

    setoutputFilepH = (value) => {
      this.setState({
        outputFilepH: value
      });
    }
    setpHmin = (value) => {
      this.setState({
        pHmin: value
      });
      console.log(this.state.pHmin)
    }
    setpHmax = (value) => {
      this.setState({
        pHmax: value
      });
    }
    setpHstep = (value) => {
      this.setState({
        pHstep: value
      });
    }
    setproteinDieletric = (value) => {
      this.setState({
        proteinDielectric: value
      });
    }
    setsolventDielectric = (value) => {
      this.setState({
        solventDielectric: value
      });
    }
    setionicStrength = (value) => {
      this.setState({
        ionicStrength: value
      })
    }

    retrieveNumberOfSites = (event) => {
      const inputFile = event.target.files[0]   
      
      console.log(inputFile)


      const numberOfSites = 1
      this.setState({
        inputFile: inputFile,
        numberOfTitratableSites: numberOfSites
    });
    }
    
    retrieveFromPDB = (object, pdbid) => {
      
      if (!pdbid) {
        object.setState({
          pdbcode_error: null
        });
        return
      }

      checkPdbId(object, pdbid)
    }
    
    submitForm = () => {
      console.log(this.state);
      var triggerHappy = true

      // Check if pdbfile has been saved
      if (!this.state.pdbfile) {
        triggerHappy = false
        console.log('pdbfile not found')
      }

      // Check that if the input mode is a uploaded file
      // we have saved the naming scheme
      if (this.state.inputModeRadio === 'upload') {
        if (!this.state.inputFileNamingScheme) {
            triggerHappy = false
            console.log('inputFileNamingScheme not found')
        }
      }

      // Check there is at least one of the output modes
      if (!(this.state.outputpKValues || this.state.outputPDBFile)) {
        triggerHappy = false
        console.log('no output mode found')
      }

      // Check that if the output file mode is turned on
      // its desired pH value and naming schema are selected
      if (this.state.outputPDBFile) {
        if (!(this.state.outputFilepH && this.state.outputFileNamingScheme)) {
          triggerHappy = false
          console.log('outputfile pH value or naming schema not found')
        }
      }

      if (this.state.outputpKValues) {
        if (isNaN(this.state.pHmin) || isNaN(this.state.pHmax) || isNaN(this.state.pHstep)) {
          triggerHappy = false
          console.log('pH values not defined', isNaN(this.state.pHmin), isNaN(this.state.pHmax), isNaN(this.state.pHstep))
        }
      }

      if (!(this.state.proteinDielectric && this.state.solventDielectric && this.state.ionicStrength)) {
        triggerHappy = false
        console.log('either dielectrics or ionic strength is not defined')
      }

      if (!triggerHappy) {
        toast.warning('Please address the missing parameters highlighted in yellow', {
          position: "top-right",
          autoClose: 5000,
          hideProgressBar: true,
          closeOnClick: true,
          pauseOnHover: true,
          draggable: true
        });
      } else {
        submit_pypka_calculation(this)
      }
    }

    render() {

      return (
    <Layout>
   
      <SEO title="PypKa Calculations Server" />
      <HeaderBar image={comp_cluster} title={"PypKa Server"} subtitle={"Leverage the webserver to quickly start calculating p<em>K</em><sub>a</sub> values in your protein of interest"} />
  
      
      <section className="section bg-gray" id="section-form">
          <div className="container">
            <header className="section-header">
              <small>input</small>
              <h2>{this.props.isDarKMode} Parameters</h2>
            </header>
  
            <div className="row gap-y">
  
              <div className="col-12 col-lg-4">
                  <div className="form-group">
                  <h3 className="divider" style={{ fontSize: "1rem", width: "80%", color: "#777"}}>Input Mode</h3>
                      <div className="custom-controls-stacked">
                          <label className="custom-control custom-radio" 
                                 style={{ textAlign: "center", margin: "auto"}}>
                                <input type="radio" id="input-mode-pdbcode"
                                       className="custom-control-input" name="radio_input_mode" val="pdb_code" 
                                       checked={this.state.inputModeRadio === "pdb_code" }
                                       onChange={(e) => {
                                          document.getElementsByName('pdb_code')[0].value = ''
                                          this.setState({
                                            pdbfile: null,
                                            inputModeRadio: e.target.getAttribute('val'),
                                            inputFileNamingScheme: "GROMOS"
                                           })
                                       }} />
                              <span className="custom-control-indicator"></span>
                              <span className="custom-control-description">PDB Code</span>
                          </label>
                          <label className="custom-control custom-radio" 
                                 style={{ textAlign: "center", margin: "auto"}}>
                                <input type="radio" id="input-mode-upload"
                                       className="custom-control-input" name="radio_input_mode" val="upload" 
                                       checked={this.state.inputModeRadio === "upload"}
                                       onChange={(e) => {
                                        document.getElementsByName('file_uploaded')[0].value = ''
                                        this.setState({
                                          pdbfile: null,
                                          inputModeRadio: e.target.getAttribute('val')
                                        })
                                    }} />
                            <span className="custom-control-indicator"></span>
                            <span className="custom-control-description">Upload File</span>
                          </label>
                      </div>
                  </div>
              </div>
  
              <div className="col-12 col-lg-4">
                  <div className="form-group"  style={{ textAlign: "center"}}>
                  <h3 className="divider" 
                      style={{ fontSize: "1rem", width: "80%", color: "#777"}}>PDB Code</h3>
                      <input className={(this.state.inputModeRadio === 'pdb_code' && [null, ''].includes(this.state.pdbcode)) ? styles.borderPrimary + ' ' + styles.formControl : styles.formControl}
                             type="text" placeholder="Ex: 4LZT" name="pdb_code" 
                             onChange={(e) => {
                                 const pdbid = e.target.value
                                 this.retrieveFromPDB(this, pdbid)
                                }}
                             disabled={this.state.inputModeRadio !== 'pdb_code'}/>
                      
                      <div className={(this.state.inputModeRadio === 'pdb_code' && ![null, ''].includes(this.state.pdbcode_error)) ? "alert alert-warning" : ""}>
                        {this.state.pdbcode_error}
                      </div>
                  </div>
                  <div className="form-group"  style={{ textAlign: "center"}}>
                  <h3 className="divider" 
                      style={{ fontSize: "1rem", width: "80%", color: "#777"}}>Input File</h3>
                      <input className={(this.state.inputModeRadio === 'upload' && [null, ''].includes(this.state.pdbfile)) ? styles.borderPrimary + ' ' + styles.formControl : styles.formControl}
                             type="file" name="file_uploaded" disabled={this.state.inputModeRadio !== 'upload'}                       
                             onChange={
                              (e) => {
                                var self = this
                                var reader = new FileReader()
                                reader.onload = function () {
                                  const content = reader.result
                                  
                                  async function reportTitrableSites(content) {
                                    console.log(content)
                                    const data = await getNumberOfTitrableSites(content)
                                    console.log(data)

                                    console.log('nchains' in data)
                                    if ('nchains' in data) {
                                      var nchains = data.nchains
                                      var nsites = data.nsites
                                    } else {
                                      var nchains = 0
                                      var nsites = 0
                                      content = ''
                                    }

                                    self.setState({
                                      pdbfile: content,
                                      nchains: nchains,
                                      nsites: nsites                                    
                                    })

                                  }
                                  reportTitrableSites(content)
                                };
                                var input_file = e.target.files[0]
                                var sys_name = input_file.name.split('.')[0]
                                this.setState({
                                  protein_name: sys_name
                                })
                                
                                reader.readAsBinaryString(input_file);
                             }} />
                      <div style={{ marginBottom: "5px" }}></div>
                      <select className={(this.state.inputModeRadio === 'upload' && [null, ''].includes(this.state.inputFileNamingScheme)) ? styles.borderPrimary + ' ' + styles.formControl : styles.formControl}
                              id="file_uploaded_force_field" 
                              disabled={this.state.inputModeRadio !== 'upload'} defaultValue={""}
                              onBlur={(e) => {
                                      console.log(e.target.value)
                                      this.setState({
                                          inputFileNamingScheme: e.target.value
                                      })}}>
                          <option value="" disabled>Select the naming scheme</option>
                          <option value="AMBER">AMBER</option>
                          <option value="CHARMM">CHARMM</option>
                          <option value="GROMOS">GROMOS</option>
                          <option value="GROMOS">PDB</option>
                      </select>
                  </div>
              </div>
  
              <div className="col-12 col-lg-4" style={{ borderLeft: "0.5px solid rgba(47, 47, 47, 0.8)" }}>
                <div className="form-group">
                  <label>Help</label>
                  <p>
                      The input PDB file can be retrieved from the Protein Data Bank or it can be uploaded by the user.
                  </p>
                  <p>
                    To use a structure from the Protein Data Bank please provide a valid 4-character PBD ID.
                  </p>
                  <p>
                    To use a uploaded structure please select a PDB file with one of the following naming schemes: AMBER, GROMOS, CHARMM or PDB.
                  </p>
                </div>
              </div>
            </div>
  
            <br /><br /><br />
            <div className="row gap-y">
  
              <div className="col-12 col-lg-4">
                  <div className="form-group">
                  <h3 className="divider" 
                      style={{ fontSize: "1rem", width: "80%", color: "#777"}}>Ouput</h3>
                      <div className="custom-controls-stacked">
                          <label className="custom-control custom-checkbox" style={{margin: "auto"}}>
                            <input className="custom-control-input"
                                   type="checkbox" checked={this.state.outputpKValues}
                                   onChange={(e) => {
                                     console.log(this.state.outputpKValues + this.state.outputPDBFile)
                                        this.setState({
                                            outputpKValues: e.target.checked
                                        })}} />
                            <span className="custom-control-indicator"></span>
                            <span className="custom-control-description" style={{width: "100px"}} >pK<sub>a</sub> values</span>
                          </label>
                          <label className="custom-control custom-checkbox" style={{margin: "auto"}}>
                              <input type="checkbox" className="custom-control-input" 
                                     onChange={(e) => {
                                        this.setState({
                                            outputPDBFile: e.target.checked
                                        })
                                        
                                      }} />
                              <span className="custom-control-indicator"></span>
                              <span className="custom-control-description" style={{width: "100px"}} >MD-ready PDB</span>
                            </label>
  
                        </div>
                        <div className={(this.state.outputpKValues + this.state.outputPDBFile === 0) ? "alert alert-warning" : styles.hidden}>
                          One of the two must be selected
                        </div>
                        
                  </div>
              </div>
  
              <div className="col-12 col-lg-4">
                  <div className="form-group" style={{ textAlign: "center"}}>
                  <h3 className="divider" 
                      style={{ fontSize: "1rem", width: "80%", color: "#777"}}>Output PDB</h3>
                      
                      <select className={(this.state.outputPDBFile === true && [null, ''].includes(this.state.outputFileNamingScheme)) ? styles.borderPrimary + ' ' + styles.formControl : styles.formControl}
                              disabled={!this.state.outputPDBFile} 
                              defaultValue={""} onBlur={(e) => {
                                this.setState({
                                  outputFileNamingScheme: e.target.value
                                })
                              }}><option value="" disabled>Select the desired naming scheme</option>
                          <option value="amber">AMBER</option>
                          <option value="gromos_cph">GROMOS</option>
                      </select>
                      <br />

                      <div className={styles.formControl}>
                        <InputFloat 
                            name={"pH value"}
                            extraStyle={"form-control"}
                            defaultValue={this.state.outputFilepH}
                            condition1={this.state.outputPDBFile === true}
                            condition2={[null, '', NaN].includes(this.state.outputFilepH)}
                            disableOn={this.state.outputPDBFile}
                            saveState={this.setoutputFilepH}
                          />
                      </div>
                  </div>
              </div>
  
              <div className="col-12 col-lg-4" style={{ borderLeft: "0.5px solid rgba(47, 47, 47, 0.8)" }}>
                <div className="form-group">
                  <p>
                      The available outputs are the p<em>K</em><sub>a</sub> values of the titratable amino acids and a PDB formatted file with the most likely protonation state at a given pH value.
                  </p>
                  <p>
                    Please select the desired output file naming scheme and the target pH value.
                  </p>
                </div>
              </div>
            </div>
  
            <br /><br /><br />
            <div className="row gap-y">
              <div className="col-8" style={{padding: "0px", marginBottom: "-15px"}}>
                <h3 className="divider" 
                    style={{ fontSize: "1rem", width: "80%", color: "#777"}}>Parameters</h3>
              </div>
              <div className="col-4">

              </div>
  
              <div className="col-12 col-lg-4">
                  <div className="form-group">

                    <div className="row">
                      <div className="col-12 col-lg-6">
                            <label>pH Titration Range</label>
                      </div>    
                      <div className="col-12 col-lg-6">
                        <div className={styles.inputGroup}>
                            <InputFloat 
                              name={"pH min"}
                              extraStyle={"form-control"}
                              defaultValue={this.state.pHmin}
                              condition1={true}
                              condition2={[null, '', NaN].includes(this.state.pHmin)}
                              disableOn={this.state.outputpKValues}
                              saveState={this.setpHmin}
                            />
                        </div>
                        
                        <div className={styles.inputGroup}>
                            <InputFloat
                              name={"pH max"}
                              extraStyle={"form-control"}
                              defaultValue={this.state.pHmax}
                              condition1={true}
                              condition2={[null, '', NaN].includes(this.state.pHmax)}
                              disableOn={this.state.outputpKValues}
                              saveState={this.setpHmax}
                            />
                        </div>
                        
                        <div className={styles.inputGroup}>
                            <InputFloat 
                              name={"pH step"}
                              extraStyle={"form-control"}
                              defaultValue={this.state.pHstep}
                              condition1={true}
                              condition2={[null, '', NaN].includes(this.state.pHstep)}
                              disableOn={this.state.outputpKValues}
                              saveState={this.setpHstep}
                            />
                        </div>
                      </div>  
                    </div>
                  </div>
              </div>
  
              <div className="col-12 col-lg-4">
                  <div className="row">
                    <div className="col-12 col-lg-6">
                        <label>Dielectric Constants</label>
                    </div>
                    <div className="col-12 col-lg-6">
                        <div className={styles.inputGroup}>
                            <InputFloat 
                                  name={"Protein"}
                                  extraStyle={"form-control"}
                                  defaultValue={this.state.proteinDielectric}
                                  condition1={true}
                                  condition2={[null, '', NaN].includes(this.state.proteinDielectric)}
                                  disableOn={true}
                                  saveState={this.setproteinDieletric}
                                />
                        </div>

                        <div className={styles.inputGroup} style={{ marginTop: "5px"}}>
                                <InputFloat 
                                  name={"Solvent"}
                                  extraStyle={"form-control"}
                                  defaultValue={this.state.solventDielectric}
                                  condition1={true}
                                  condition2={[null, '', NaN].includes(this.state.solventDielectric)}
                                  disableOn={true}
                                  saveState={this.setsolventDieletric}
                                />
                        </div>
                    </div>
                  </div>
                  <div className="row"  style={{ marginTop: "15px"}}>
                    <div className="col-12 col-lg-6">
                      <label>Ionic Strength</label>
                    </div>
                    <div className="col-12 col-lg-6">
                      <InputFloatSimple 
                        width={"100%"}
                        extraStyle={"form-control"}
                        defaultValue={this.state.ionicStrength}
                        condition1={true}
                        condition2={[null, '', NaN].includes(this.state.ionicStrength)}
                        disableOn={true}
                        saveState={this.setionicStrength}
                      />

                    </div>
                  </div>
              </div>
  
              <div className="col-6 col-lg-4" style={{borderLeft: "0.5px solid rgba(47, 47, 47, 0.8)" }}>
                <div className="form-group">
                  <p>
                    Please set pH min and pH max so that the pKa of all titratable amino acids is withing the titration range. 
                    pH step determines the difference between consecutive pH values in the simulated titration.
                  </p>
                  <p>
                    The default value for the protein dielectric has been found to minimize the error of the prediction. 
                    The default value for the solvent dielectric is the experimental value of water at around 273K.
                  </p>
                </div>
              </div>
            </div>
          </div>
      </section>
      <section className="section section-inverse py-40" style={{ backgroundColor: "#8ea6e6" }}>
          <div className="container">
            <div className="row gap-y align-items-center text-center text-md-left" style={{ alignItems: "flex-end!important" }}>
              <div className="col-12 col-md-9">
                <h4 className="fw-300 mb-0">{this.state.nsites} Titratable Sites over {this.state.nchains} Chains</h4>
                            <h4 className="fw-300 mb-0">Estimated Running Time: {getTimeEstimate(this.state.nsites, this.state.pHmax, this.state.pHmin, this.state.pHstep)}s</h4>
                <h4 className={styles.hidden + " fw-300 mb-0"} >Queued Submissions: 0 (Average Time per Submission: 60s)</h4>
              </div>
  
              <div className="col-12 col-md-3 text-center text-md-right">
                <a className="btn btn-lg btn-round btn-white" onClick={this.submitForm} style={{ color: "#b5b9bf"}}> Submit</a>
              </div>
            </div>
          </div>
      </section>
      
    </Layout>
)}
        }

export default RunPage