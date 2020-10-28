from flask import Flask, render_template, request, redirect
from werkzeug.utils import secure_filename


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
        hssp = request.form['hssp']
        query = request.form['query_seq']
        w_l = request.form['word_length']
        insertion = request.form['insertion']
        deletion = request.form['deletion']
        mismatch = request.form['mismatch']
        mscore = request.form['mscore']
        print(hssp, query, w_l, insertion, deletion, mismatch, mscore)

        return render_template('results.html')

if __name__ == '__main__':
    app.run()