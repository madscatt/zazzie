import json
import os
import sys


def stuhrmann_plot(json_filename):

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
1/contrast, Rg^2, Rg^2 error, calculated Rg^2, Rg^2 - calculated Rg^2 [delta_rho_inverse, rg_squared, rg_squared_error, rg_squared_calculated, diff
        """

    lineplot = {
        "data": [
            {
                "x": json_object["1/contrast"],
                "y": json_object["calculated_Rg^2"],
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
                "hoverinfo": "x+y",
            },
            {
                "x": json_object["1/contrast"],
                "y": json_object["Rg^2"],
                "xaxis": "x",
                "yaxis": "y",
                "error_y": {
                    "type": "data",
                    "array": json_object["Rg^2_error"],
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
                "x": json_object["1/contrast"],
                "y": json_object["Rg^2-calculated_Rg^2"],
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
            }
        ],
        "layout": {
            "height": 700,
            "grid": {
                "rows": 2,
                "columns": 1,
                "subplots": [["xy"], ["x2y2"]],
            },
            "title": {
                "text": "Stuhrmann plot",
                "y": 0.9,
                "x": 0.47,
            },
            "hovermode": "x",
            "xaxis": {
                "title": {
                    "text": "1/contrast",
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
                    "text": "1/contrast",
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
                "title": "R<sub>g</sub> <sup>2</sup>",
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
    json_filename = "general_output_file_stuhrmann.out.json"
    lineplot = stuhrmann_plot(json_filename)
#    print(lineplot)
    open("stuhrmann_plot_for_genapp.json", "w").write(
        json.dumps(lineplot, indent=2))
