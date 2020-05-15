import { Link } from "gatsby"
import React from "react"

const Footer = () => (
    <footer className="site-footer">
      <div className="container">
        <div className="row gap-y">
          <div className="col-12 col-md-6">
            <p className="text-center text-md-left">PypKa 2020 <a className="text-dark" href="http://mms.rd.ciencias.ulisboa.pt">MMS@FCUL</a>.
              Crafted by <a className="text-dark" href="https://github.com/pedrishi">Pedro Reis</a>
            </p>
          </div>

          <div className="col-12 col-md-6">
            <ul className="nav nav-inline nav-primary nav-dotted nav-dot-separated justify-content-center justify-content-md-end">
              <li className="nav-item">
                <Link className="nav-link" to="/">Terms of use</Link>
              </li>
              <li className="nav-item">
                <Link className="nav-link" to="/">Privacy policy</Link>
              </li>
            </ul>
          </div>
        </div>
      </div>
    </footer>
)


export default Footer