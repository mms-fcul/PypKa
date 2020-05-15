import { Link } from "gatsby"
import React from "react"
import logo from "../images/pypka_logo.png"
import logo_inv from "../images/pypka_logo_inv.png"
import { node } from "prop-types"

const Header = ({ siteTitle }) => (
    <nav className="topbar topbar-inverse topbar-expand-md topbar-sticky">
        <div className="container">
          <div className="topbar-left">
            <button className="topbar-toggler">&#9776;</button>
            <Link className="topbar-brand" to="/">
              <img className="logo-default" src={ logo } alt="logo" />
              <img className="logo-inverse" src={ logo_inv } alt="logo" />
            </Link>
          </div>
          <div className="topbar-right">
            <ul className="topbar-nav nav">
              <li className="nav-item"><Link className="nav-link" to="/">Home</Link></li>
              <li className="nav-item"><Link className="nav-link" to="/run-pypka">Run Pypka</Link></li>
              <li className="nav-item" style={{display: "none"}}><Link className="nav-link" to="/latest">Latest Simulations</Link></li>
              <li className="nav-item"><a className="nav-link" href="https://github.com/mms-fcul/PypKa">Code</a></li>
              <li className="nav-item"><a className="nav-link" href="https://pypka.readthedocs.io/en/latest/">Documentation</a></li>
              <li className="nav-item"><a className="nav-link" href="http://mms.rd.ciencias.ulisboa.pt">MMS@FCUL</a></li>
            </ul>
          </div>
        </div>
    </nav>
)


export default Header
