import json
import os
import sys


def decomposition_plots(json_filename):

    print("json_filename: ", json_filename)
    if os.path.isfile(json_filename):
        print("File exists: ", json_filename)
        # Read the JSON data from the file
        with open(json_filename, 'r') as f:
            json_object = json.load(f)
#            print("json_object: ", json_object)
    else:
        error = "File does not exist: " + json_filename
        print(error)
        sys.exit()

# data is a list of dictionaries, each representing a plot for one item in 'scattering_data' or 'rescaled_scattering_data'
    data = []
    for i in range(len(json_object['fraction_d2o'])):
        data.append(
            {
                "x": json_object['q_values'],
                "y": json_object['scattering_data'][i],
                "xaxis": "x",
                "yaxis": "y",
                "error_y": {
                    "type": "data",
                    "array": json_object['scattering_data_error'][i],
                    "visible": True,
                    "thickness": 1.5,
                    "width": 0,
                },
                "type": "scatter",
                "mode": "markers",
                "name": str(json_object['fraction_d2o'][i]),
                "marker": {
                    "size": 4.5,
                    "symbol": "circle"
                },
                "showlegend": True,
                "hoverinfo": "x+y",
            }
        )
        data.append(
            {
                "x": json_object['q_values'],
                "y": json_object['rescaled_scattering_data'][i],
                "xaxis": "x",
                "yaxis": "y",
                "error_y": {
                    "type": "data",
                    "array": json_object['rescaled_scattering_data_error'][i],
                    "visible": True,
                    "thickness": 1.5,
                    "width": 0,
                },
                "type": "scatter",
                "mode": "markers",
                "name": str(json_object['fraction_d2o'][i]) + " rescaled",
                "marker": {
                    "size": 4.5,
                    "symbol": "x"
                },
                "showlegend": True,
                "hoverinfo": "x+y",
            }
        )
#    print("data: ", data)

    lineplot1 = {
        "data": data,
        "layout": {
            #            "height": 700,
            "legend": {"title": {"text": 'Fraction D<sub>2</sub>O', }, },
            "colorway": ['red', 'red', 'blue', 'blue', 'black', 'black', 'green', 'green', 'magenta', 'magenta', 'cyan', 'cyan', 'lime', 'lime', 'purple', 'purple', 'orange', 'orange', 'grey', 'grey'],
            "title": {
                "text": "Experimental and Rescaled Scattering Data",
                "y": 0.9,
                "x": 0.47,
            },
            "hovermode": "x",
            "xaxis": {
                "type": "log",
                "title": {
                    "text": "q (Å<sup>-1</sup>)",
                },
                "autorange": True,
                "showgrid": True,
                "zeroline": False,
                "showline": True,
                "autotick": True,
                "ticks": "inside",
                "mirror": "ticks",
                "showticklabels": True,
                "hoverformat": ".2r",
            },
            "yaxis": {
                "type": "log",
                "title": "I(q)",
                "autorange": True,
                "showgrid": True,
                "zeroline": False,
                "showline": True,
                "autotick": True,
                "ticks": "inside",
                "mirror": "ticks",
                "showticklabels": True,
                "exponentformat": "e",
                "hoverformat": ".2r",
            },
        },

        "config": {
            "responsive": True,
            "scrollZoom": True,
            "toImageButtonOptions": {
                "format": "png",
                "scale": 2,
            },
            "modeBarButtonsToRemove": ["select2d", "lasso2d"],
            "modeBarButtonsToAdd": ["togglespikelines", "v1hovermode", "hovercompare"],
            "showLink": True,
            "plotlyServerURL": "https://chart-studio.plotly.com",
        }
    }

    return lineplot1


if __name__ == "__main__":
    json_filename = "general_output_file.out.json"
    lineplot1 = decomposition_plots(json_filename)
#    print(lineplot1)
# this JSON file can be pasted into the GenApp advanced tutorial plotly test utility
    open("decomposition_plot1_for_genapp.json", "w").write(
        json.dumps(lineplot1, indent=2))
