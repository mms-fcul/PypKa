import React from "react"
import { Link } from "gatsby"

import Layout from "../components/layout"
import SEO from "../components/seo"
import library from "../images/library.jpg"
import { graphql } from 'gatsby'
import CustomParticles from "../components/particles.js";

const IndexPage = ({ data }) => {
return (
  <Layout>
    <SEO title="Home" />
    
    <header className="header header-inverse h-fullscreen pb-80" style={{backgroundColor: "#4f407b"}}>

    <CustomParticles />
      <div className="container text-center">

        <div className="row h-full">
          <div className="col-12 col-lg-8 offset-lg-2 align-self-center">

            <h1 className="display-1 hidden-sm-down">PypKa Server</h1>
            <h1 className="display-4 hidden-md-up">PypKa Server</h1>
            <br />
            <p className="fs-20 w-400 mx-auto hidden-sm-down">Flexible Poisson-Boltzmann based pK<sub>a</sub> calculations with proton tautomerism</p>
            <p className="fs-16 w-250 mx-auto hidden-md-up">Flexible Poisson-Boltzmann based pK<sub>a</sub> calculations with proton tautomerism</p>

            <hr className="w-80" />
            <br />

            <Link className="btn btn-xl btn-round btn-white w-200 hidden-sm-down" to="/run-pypka">Run PypKa</Link>
            <Link className="btn btn-lg btn-round btn-white w-200 hidden-md-up" to="/run-pypka">Run PypKa</Link>

          </div>

          <div className="col-12 align-self-end text-center">
            <Link className="scroll-down-1 scroll-down-inverse" to="/#why" data-scrollto="section-intro"><span></span></Link>
          </div>

        </div>

      </div>
    </header>

    <main className="main-content">
    </main>


    <section className="section" style={{margin: "30px 50px 100px",
    boxShadow: "0 0 35px rgba(0, 0, 0, 0.1)",
    border: "1px solid #ebebeb"}}
    >
      <div className="container">
        <header className="section-header">
          <small>Features</small>
          <h2>Why should you use PypKa?</h2>
          <hr />
          <p className="lead"></p>
        </header>

        <div className="row gap-y">

          <div className="col-12 col-md-6 col-xl-4 feature-1">
            <p className="feature-icon text-dark"><i className="icon-grid"></i></p>
            <h5>Poisson–Boltzmann-based</h5>
            <p>Electrostatic-driven p<em>K</em><sub>a</sub> estimations</p>
          </div>


          <div className="col-12 col-md-6 col-xl-4 feature-1">
            <p className="feature-icon text-info"><i className="icon-search"></i></p>
            <h5>Open-source</h5>
            <p>Inspect the source code and hack it as you see fit.</p>
          </div>


          <div className="col-12 col-md-6 col-xl-4 feature-1">
            <p className="feature-icon"><i className="icon-attachment"></i></p>
            <h5>DNA, Membrane & Ions</h5>
            <p>Include background charges in your calculations.</p>
          </div>


          <div className="col-12 col-md-6 col-xl-4 feature-1">
            <p className="feature-icon text-warning"><i className="icon-layers"></i></p>
            <h5>Multiprocessing Capabilities</h5>
            <p>Exploit the linear scaling of pypka by deploying it in your computational cluster.</p>
          </div>


          <div className="col-12 col-md-6 col-xl-4 feature-1">
            <p className="feature-icon text-danger"><i className="icon-recycle"></i></p>
            <h5>Reusable & extensible</h5>
            <p>Integrate it with your own software or scripts.</p>
          </div>


          <div className="col-12 col-md-6 col-xl-4 feature-1">
            <p className="feature-icon text-success"><i className="icon-adjustments"></i></p>
            <h5>Flexible</h5>
            <p>All parameters used in the calculations are tweackable in the input data.</p>
          </div>

        </div>

      </div>
    </section>

    <section className="section py-150" style={{backgroundImage: `url(${library})`, backgroundSize: "cover"}} data-overlay="7" id="why">
      <div className="container">

        <div className="section-dialog text-center">
          <h2>Citation</h2>
          <p>
            Are you using PypKa in your research? Please cite:<br />
            <code>Pypka article in JSIM</code>
          </p>

          </div>

      </div>
    </section>

    <section className="section pb-0 overflow-hidden hidden-sm-down">
      <div className="container">
        <header className="section-header">
          <small><strong>Statistics</strong></small>
          <h2>Fast and Accurate</h2>
          <hr />
          <p className="lead">Easily calculate p<em>K</em><sub>a</sub> values in your protein of interest with our server.</p>
        </header>

        <div className="row gap-y text-center">

          <div className="col-md-4 d-flex flex-column">
            <div className="mb-60">
              <p className="text-info fs-50 mb-0">2</p>
              <p>Force Field Choices</p>
            </div>
            <div className="px-20 mt-auto">
              <img className="shadow-4 opacity-80 aos-init aos-animate" src={ data.tabfiles.edges[0].node.publicURL } alt="..." data-aos="slide-up" data-aos-delay="300" />
            </div>
          </div>

          <div className="col-md-4 d-flex flex-column">
            <div className="mb-7">
              <span className="text-info fs-50">40s</span><br />
              <p>Calculation on a 50 residue Protein</p>
            </div>
            <div className="mt-auto">
              <img className="shadow-6 aos-init aos-animate" src={ data.tabfiles.edges[1].node.publicURL } alt="..." data-aos="slide-up" />
            </div>
          </div>

          <div className="col-md-4 d-flex flex-column">
            <div className="mb-7">
              <span className="text-info fs-50">0.3</span><br />
              <p>Average Error in a 500 residue dataset</p>
            </div>
            <div className="px-20 mt-auto">
              <img className="shadow-4 opacity-80 aos-init aos-animate" src={ data.tabfiles.edges[2].node.publicURL } alt="..." data-aos="slide-up" data-aos-delay="600" />
            </div>
          </div>

        </div>

      </div>
    </section>

    <section className="section bg-gray text-center">

      <div className="row">
        <div className="col-12 offset-md-3 col-md-6">
          <p><img src={ data.githubimg.edges[0].node.publicURL } alt="..." /></p>
          <br />
          <h3 className="fw-900 mb-20">Fork the project on GitHub</h3>
          <p className="lead text-muted">Pypka is open source! It's developed and maintained on GitHub so that you can make it your own.</p>
          <br />
          <a className="btn btn-lg btn-round btn-success" href="https://github.com/mms-fcul/PypKa">View GitHub Project</a>
        </div>
      </div>

    </section>
    
  </Layout>
)}

export const query = graphql`
{
  tabfiles: allFile(filter: {relativePath: {regex: "/tab/"}}) {
    edges {
      node {
        publicURL
      }
    }
  }
	githubimg:allFile(filter: {relativePath: {regex: "/git/"}}) {
    edges {
      node {
        publicURL
      }
    }
  }  
}
`


export default IndexPage
