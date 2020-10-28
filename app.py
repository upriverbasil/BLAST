from flask import Flask, render_template, request, redirect
from werkzeug.utils import secure_filename
import BLAST


app = Flask(__name__)


@app.route('/')
def index():
    return render_template("index.html")

@app.route('/upload', methods=['POST'])
def upload_file():
    uploaded_file = request.files['file']
    if uploaded_file.filename != '':
        uploaded_file.save("db_seq.fasta")
    print(uploaded_file.filename)
    return redirect(request.referrer)


@app.route('/submit', methods=['POST'])
def submit():
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