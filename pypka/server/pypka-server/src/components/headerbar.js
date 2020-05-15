import React from "react"


const HeaderBar = (props) => (
<header className="header header-inverse bg-fixed" style={{ backgroundImage: `url(${props.image})` }} data-overlay="8">
    <div className="container text-center">

      <div className="row">
        <div className="col-12 col-lg-8 offset-lg-2">

          <h1>{props.title}</h1>
          <p className="fs-18 opacity-70" dangerouslySetInnerHTML={{__html: props.subtitle}}></p>

        </div>
      </div>

    </div>
</header>
)

export default HeaderBar