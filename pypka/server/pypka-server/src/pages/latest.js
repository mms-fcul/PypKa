import React from "react"
import { Link } from "gatsby"

import comp_cluster from "../images/comp_cluster.jpg"

import Layout from "../components/layout"
import HeaderBar from "../components/headerbar"
import ControlledExpansionPanels from "../components/accordion"

import styles from "./run-pypka.module.css"

import axios from "axios"

var config = { headers: {  
  'Content-Type': 'application/json',
  'Access-Control-Allow-Origin': '*'}
}

async function updateLatestSimulations() {
  try {
    const response = await axios.post('http://127.0.0.1:5000/getLatestsSubmissions', {}, config)
    const data = response.data
    console.log(data)
    return data
  } catch (error) {
    console.log(error)
  }
}

const Latest = () => {

  const componentWillMount = () => {
    updateLatestSimulations()
  }

  return (
  <Layout>
    <HeaderBar image={comp_cluster} 
             title={"Latest Simulations"} 
             subtitle={""} />

    <section id="basic" className="section bb-1">
      <div className="container">
        <div className="row gap-y">
          <div className="col-md-8 offset-md-2">
            
              <ControlledExpansionPanels />

            
          </div>
        </div>
      </div>
    </section>
  </Layout>
)}

export default Latest