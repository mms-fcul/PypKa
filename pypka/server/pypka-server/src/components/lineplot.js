import React from 'react';
import Plot from 'react-plotly.js';

class LinePlot extends React.Component {
    constructor(props) {
      super(props);
      this.state = { 
          x: props.x, 
          y: props.y,
          width: props.width
      };
    }
    
  render() {
    console.log(this.state)
    return (
      <Plot
        data={[
          {
            type: 'scatter', 
            x: this.state.x, 
            y: this.state.y,    
            mode: 'lines+markers',
            marker: {color: 'red'},
            line: {shape: 'spline'}  
          }
        ]}
        layout={ {
            title: 'Titration Curve', 
            autosize: true, 
            xaxis: {
                title: 'pH',
                showgrid: false,
                zeroline: true
            },
            yaxis: {
                title: 'Protonation',
                showline: false
            }
        } }
        style={ {width: "100%", height: "100%"} }
        config={ {displayModeBar: false, responsive: true} }        
      />
    );
  }
}

export default LinePlot