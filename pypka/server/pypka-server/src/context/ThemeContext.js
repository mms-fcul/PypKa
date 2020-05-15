class GlobalState {
    constructor(subID) {
        this.state = {
            subID: subID,
            pKas: [],
            tit_x: [],
            tit_y: [],
            params: '',
            pdb_out: null
        }

        const savedSubmission = JSON.parse(localStorage.getItem(subID))
        if (savedSubmission) {
            this.state.pKas = savedSubmission.pKas
            this.state.tit_x = savedSubmission.titration[0]
            this.state.tit_y = savedSubmission.titration[1]
            this.state.params = savedSubmission.parameters
            this.state.pdb_out = savedSubmission.pdb_out
        }
    }

    saveSubmission = (subID, results) => {
      localStorage.setItem(subID, JSON.stringify(results))
    }
}

export default GlobalState