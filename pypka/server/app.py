from flask import Flask, render_template

app = Flask(__name__)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/run-pypka')
def run_pypka():
    return render_template('runpypka.html')

@app.route('/latest-simulations')
def latest_sims():
    return render_template('latest.html')