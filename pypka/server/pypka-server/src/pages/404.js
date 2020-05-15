import React from "react"

import Layout from "../components/layout"
import HeaderBar from "../components/headerbar"
import SEO from "../components/seo"

import broken_pc from "../images/broken-pc.jpg"

const NotFoundPage = () => (
  <Layout>
    <HeaderBar image={broken_pc} title={"NOT FOUND"} subtitle={"You just hit a route that doesn&#39;t exist... the sadness."} />
  </Layout>
)

export default NotFoundPage
