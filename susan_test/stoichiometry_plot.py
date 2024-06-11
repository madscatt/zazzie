import json
import os
import sys


def stoichiometry_plot(json_filename):

    print("json_filename: ", json_filename)
    if os.path.isfile(json_filename):
        print("File exists: ", json_filename)
        # Read the JSON data from the file
        with open(json_filename, 'r') as f:
            json_object = json.load(f)
    else:
        error = "File does not exist: " + json_filename
        print(error)
        sys.exit()

        """
measured I(0), I(0) error, back-calculated I(0) using the calculated Mw values, I(0) - calculated I(0) (at each fraction D2O) [fraction_d2o, izero, izero_error, izero_calc, diff]
        """

    lineplot = {
        "data": [
            {
                "x": json_object["fraction_d2o"],
                "y": json_object["izero_calc"],
                "xaxis": "x",
                "yaxis": "y",
                "type": "scatter",
                "mode": "markers",
                "name": "calculated",
                "marker": {
                    "color": "red",
                    "size": 10,
                    "symbol": "circle"
                },
                "hoverinfo": "x+y"
            },
            {
                "x": json_object["fraction_d2o"],
                "y": json_object["izero"],
                "xaxis": "x",
                "yaxis": "y",
                "error_y": {
                    "type": "data",
                    "array": json_object["izero_error"],
                    "visible": True,
                    "thickness": 3,
                    "width": 0,
                },
                "type": "scatter",
                "mode": "markers",
                "name": "data",
                "marker": {
                    "color": "blue",
                    "size": 8,
                    "symbol": "circle"
                },
                "hoverinfo": "x+y",
            },
            {
                "x": json_object["fraction_d2o"],
                "y": json_object["I(0)-calculated_I(0)"],
                "xaxis": "x2",
                "yaxis": "y2",
                "type": "scatter",
                "mode": "markers",
                "name": "residuals",
                "showlegend": True,
                "marker": {
                    "color": "black",
                    "size": 8,
                    "symbol": "circle"
                },
                "hoverinfo": "x+y",
            },

        ],
        "layout": {
            "height": 700,
            "grid": {
                "rows": 2,
                "columns": 1,
                "subplots": [["xy"], ["x2y2"]],
            },
            "title": {
                "text": "Stoichiometry I(0) plot",
                "y": 0.9,
                "x": 0.47,
            },
            "hovermode": "x",
            "xaxis": {
                "title": {
                    "text": "fraction D<sub>2</sub>O",
                    "standoff": 5,
                },
                "margin": {
                    "autoexpand": True,
                },
                "autorange": True,
                "showgrid": True,
                "zeroline": False,
                "showline": True,
                "autotick": True,
                "ticks": "inside",
                "mirror": "ticks",
                "showticklabels": True,
                "tickformat": ".2r",
                "hoverformat": ".2r",
                "domain": [0, 1],
            },
            "xaxis2": {
                "title": {
                    "text": "fraction D<sub>2</sub>O",
                    "standoff": 5,
                },
                "margin": {
                    "autoexpand": True,
                },
                "autorange": True,
                "showgrid": True,
                "zeroline": False,
                "showline": True,
                "autotick": True,
                "ticks": "inside",
                "mirror": "ticks",
                "showticklabels": True,
                "tickformat": ".2r",
                "hoverformat": ".2r",
                "domain": [0, 1],
            },
            "yaxis": {
                "title": "I(0)",
                "autorange": True,
                "showgrid": True,
                "zeroline": False,
                "showline": True,
                "autotick": True,
                "ticks": "inside",
                "mirror": "ticks",
                "showticklabels": True,
                "tickformat": ".2r",
                "hoverformat": ".2r",
                "domain": [0.47, 1.0],
            },
            "yaxis2": {
                "title": "residuals",
                "autorange": True,
                "showgrid": True,
                "zeroline": True,
                "showline": True,
                "autotick": True,
                "ticks": "inside",
                "mirror": "ticks",
                "showticklabels": True,
                "tickformat": ".2r",
                "hoverformat": ".2r",
                "domain": [0, 0.38]
            }



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

    return lineplot


if __name__ == "__main__":
    json_filename = "general_output_file_stoichiometry.out.json"
    lineplot = stoichiometry_plot(json_filename)
#    print(lineplot)
    open("stoichiometry_plot_for_genapp.json", "w").write(
        json.dumps(lineplot, indent=2))
