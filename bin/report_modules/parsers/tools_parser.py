import json

import yaml
from pygments import highlight
from pygments.formatters import HtmlFormatter
from pygments.lexers import JsonLexer


def parse_tools_yaml():
    with open("software_versions.yml") as f:
        tools_dict = yaml.safe_load(f)
        formatted_tools_json = highlight_json(
            json.dumps(format_tools_dict(tools_dict), indent=4)
        )

    return tools_dict, formatted_tools_json


def highlight_json(json_string):
    lexer = JsonLexer()
    formatter = HtmlFormatter()

    return highlight(json_string, lexer, formatter)


def format_tools_dict(input_dict):
    output_list = []
    for top_level_value in input_dict.values():
        for key, value in top_level_value.items():
            if (key, value) not in output_list:
                output_list.append((key, value))

    return dict(sorted(output_list, key=lambda x: x[0]))
