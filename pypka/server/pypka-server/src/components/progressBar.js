import React from 'react';
import { lighten, makeStyles, withStyles } from '@material-ui/core/styles';
import CircularProgress from '@material-ui/core/CircularProgress';
import LinearProgress from '@material-ui/core/LinearProgress';

const ColorCircularProgress = withStyles({
  root: {
    color: '#0e97ff',
  },
})(CircularProgress);

const ColorLinearProgress = withStyles({
  colorPrimary: {
    backgroundColor: '#9ad2ff',
  },
  barColorPrimary: {
    backgroundColor: '#0e97ff',
  },
})(LinearProgress);

const useStyles = makeStyles((theme) => ({
  root: {
    flexGrow: 1,
    margin: 0,
    textAlign: "center"
  },
  margin: {
    margin: theme.spacing(1),
  },
}));

export class ProgressLinear extends React.Component {
    state = {
        update_step: 100 * 0.5 / this.props.run_time,
        value: 0
    }

    componentDidMount() {
        setInterval(() => { 
            this.setState({
                value: this.state.value + this.state.update_step
            })
        }, 500);
    }


    render() {
    return (
        <div className={useStyles.root}>
          <ColorLinearProgress className={useStyles.root} variant="determinate" value={this.state.value}/>
        </div>
    );
  }
}

export function ProgressCircle() {
    const classes = useStyles();
  
    return (
      <div className={classes.root}>
        <ColorCircularProgress size={30} thickness={4} />        
      </div>
    );
}