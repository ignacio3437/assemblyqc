from pathlib import Path

from jinja2 import Environment, FileSystemLoader


class ReportPrinter:
    def __init__(self):
        project_dir = "/".join(__file__.split("/")[0:-1])
        path = Path(f"{project_dir}/templates")

        self.file_loader = FileSystemLoader(path)
        self.env = Environment(loader=self.file_loader)

    def print(self, stats):
        template = self.env.get_template("base.html")
        return template.render(all_stats_dicts=stats)
