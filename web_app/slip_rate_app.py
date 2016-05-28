from flask import Flask
from flask import request
from flask import render_template

app = Flask(__name__)

@app.route('/', methods=['GET', 'POST'])
def plot_form():

    fields = {}

    if request.method == 'POST':
        print(request.form)

        unit_1 = ast.literal_eval(request.form['unit_1'])
        age_1 = ast.literal_eval(request.form['age_1'])
        age_err_1 = ast.literal_eval(request.form['age_err_1'])
        age_type_1 = ast.literal_eval(request.form['age_type_1'])
        age_units_1 = ast.literal_eval(request.form['age_units_1'])
        offset_1 = ast.literal_eval(request.form['offset_1'])
        offset_err_1 = ast.literal_eval(request.form['offset_err_1'])
        offset_type_1 = ast.literal_eval(request.form['offset_type_1'])
        offset_units_1 = ast.literal_eval(request.form['offset_units_1'])
        
        result = 'something'
    else:
        result = None

    return render_template("slip_rate_form.html", result=result)


if __name__ == '__main__':
    app.debug = True
    app.run()
