import json
import os
import sys


def fraction_d2o_plots(json_filename):

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

    lineplot = {

        "data": [
            {
                "x": json_object["fraction_d2o"],
                "y": json_object["sqrt[I(0)/c]_calc"],
                "xaxis": "x",
                "yaxis": "y",
                "type": "scatter",
                "mode": "lines",
                "name": "calculated",
                "line": {
                    "dash": "dash",
                    "color": "red",
                    "width": 3,
                },
                "hoverinfo": "x+y",
            },
            {
                "x": json_object["fraction_d2o"],
                "y": json_object["sqrt[I(0)/c]"],
                "xaxis": "x",
                "yaxis": "y",
                "error_y": {
                    "type": "data",
                    "array": json_object["sqrt[I(0)/c]_error"],
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
                    "symbol": "circle",
                },
                "hoverinfo": "x+y",
            },
            {
                "x": json_object["match_point"],
                "y": json_object["match_point_y_value"],
                "xaxis": "x",
                "yaxis": "y",
                "error_x": {
                    "type": "data",
                    "array": json_object["match_point_error"],
                    "visible": True,
                    "thickness": 3,
                    "width": 0,
                },
                "type": "scatter",
                "mode": "markers",
                "name": "match point",
                "marker": {
                        "color": "green",
                        "size": "8",
                        "symbol": "circle",
                },
                "hoverinfo": "x+y",
            },
            {
                "x": json_object["fraction_d2o"],
                "y": json_object["sqrt[I(0)/c]-sqrt[I(0)/c]_calc"],
                "xaxis": "x2",
                "yaxis": "y2",
                "type": "scatter",
                "mode": "markers",
                "marker": {
                    "color": "black",
                    "size": 8,
                    "symbol": "circle",
                },
                "name": "residuals",
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
                "text": "match point",
                        "y": 0.9,
                        "x": 0.47
            },
            "hovermode": "x",
            "xaxis": {
                "title": {
                    "text": "Fraction D<sub>2</sub>O",
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
            },
            "xaxis2": {
                "title": {
                    "text": "Fraction D<sub>2</sub>O",
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
            },
            "yaxis": {
                "title": "sqrt[I(0)/c]",
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
                "domain": [0, 0.38],
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
        },

    }

    return lineplot


if __name__ == "__main__":
    json_filename = "general_output_file_match_point.out.json"  
    lineplot = fraction_d2o_plots(json_filename)
#    print(lineplot)
    open("match_point_plot_for_genapp.json", "w").write(
        json.dumps(lineplot, indent=2))
