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
    else:
        error = "File does not exist: " + json_filename
        print(error)
        sys.exit()

    lineplot3 = {
        "data": [
            {
                "x": json_object['q_values'],
                "y": json_object['composite_scattering_intensity_1'],
                "xaxis": "x",
                "yaxis": "y",
                "error_y": {
                    "type": "data",
                    "array": json_object['composite_scattering_intensity_1_error'],
                    "visible": True,
                    "thickness": 1.5,
                    "width": 0,
                },
                "type": "scatter",
                "mode": "markers",
                "name": str(json_object['component_name'][0]),
                "marker": {
                    "color": "red",
                    "size": 4.5,
                    "symbol": "circle"
                },
                "showlegend": True,
                "hoverinfo": "x+y",
            },
            {
                "x": json_object['q_values'],
                "y": json_object['composite_scattering_intensity_2'],
                "xaxis": "x",
                "yaxis": "y",
                "error_y": {
                    "type": "data",
                    "array": json_object['composite_scattering_intensity_2_error'],
                    "visible": True,
                    "thickness": 1.5,
                    "width": 0,
                },
                "type": "scatter",
                "mode": "markers",

                "name": str(json_object['component_name'][1]),
                "marker": {
                    "color": "blue",
                    "size": 4.5,
                    "symbol": "circle"
                },
                "showlegend": True,
                "hoverinfo": "x+y",
            },
            {
                "x": json_object["q_values"],
                "y": json_object["composite_scattering_intensity_12"],
                "xaxis": "x2",
                "yaxis": "y2",
                "error_y": {
                    "type": "data",
                    "array": json_object['composite_scattering_intensity_12_error'],
                    "visible": True,
                    "thickness": 1.5,
                    "width": 0,
                },
                "type": "scatter",
                "mode": "markers",
                "name": str(json_object['component_name'][0] + '-' + json_object['component_name'][1] + '\n' + 'cross term'),
                "showlegend": True,
                "marker": {
                    "color": "black",
                    "size": 4.5,
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
                "text": "Composite Scattering Intensities",
                "y": 0.9,
                "x": 0.47,
            },
            "hovermode": "x",
            "xaxis": {
                "type": "log",
                "title": {
                    "text": "q (Å<sup>-1</sup>)",
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
                #                "tickformat": ".2r",
                "hoverformat": ".2r",
                "domain": [0, 1],
            },
            "xaxis2": {
                "title": {
                    "text": "q (Å<sup>-1</sup>)",
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
                "type": "log",
                "title": "I(q)",
                "autorange": True,
                "exponentformat": "e",
                "showgrid": True,
                "zeroline": False,
                "showline": True,
                "autotick": True,
                "ticks": "inside",
                "mirror": "ticks",
                "showticklabels": True,
                #                "tickformat": ".2r",
                "hoverformat": ".2r",
                "domain": [0.47, 1.0],
            },
            "yaxis2": {
                "title": "I(q)",
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

    return lineplot3


if __name__ == "__main__":
    json_filename = "general_output_file.out.json"
    lineplot3 = decomposition_plots(json_filename)
#    print(lineplot3)
# this JSON file can be pasted into the GenApp advanced tutorial plotly test utility
    open("decomposition_plot3_for_genapp.json", "w").write(
        json.dumps(lineplot3, indent=2))
