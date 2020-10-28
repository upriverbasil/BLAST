from flask import Flask, render_template, request, redirect
from werkzeug.utils import secure_filename
import BLAST

# Initiate Flask object with the current __name__ env variable
app = Flask(__name__)


@app.route('/') # Load index.html on default
def index():
    """
    Renders index.html
    """
    return render_template("index.html")

@app.route('/upload', methods=['POST']) # Executed when uploading database fasta file
def upload_file():
    """
    Store the uploaded fasta file as db_seq.fasta
    """
    uploaded_file = request.files['file']
    if uploaded_file.filename != '':
        uploaded_file.save("db_seq.fasta")
    print(uploaded_file.filename)
    return redirect(request.referrer)


@app.route('/submit', methods=['POST']) # Executed when results are submitted
def submit():
    """
    stores hyperparameters from form and calls the BLAST constructor
    runs and renders the result on result.html
    """
    if request.method == "POST":
        hssp = int(request.form['hssp'])
        query = str(request.form['query_seq'])
        w_l = int(request.form['word_length'])
        insertion = int(request.form['insertion'])
        deletion = int(request.form['deletion'])
        mismatch = int(request.form['mismatch'])
        mscore = int(request.form['mscore'])
        blast = BLAST.BLAST("db_seq.fasta", query, w_l, hssp, insertion,deletion,mismatch,mscore)
        alg = blast.run()
        return render_template('result.html', alignments = alg)

if __name__ == '__main__':
    app.run()