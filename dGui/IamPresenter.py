import IamSettingsBox
import os
import json
import subprocess

class IamPresenter:
    def __init__(self, iam_view: IamSettingsBox):
        self.iam_view = iam_view
        self.radiation = "X-ray"
        self.table = self.iam_view.current_table()
        self.radiation = self.iam_view.current_radiation()

    def save(self, directory):
        model = {
            "model": "IAM",
            "electron scattering": self.radiation == "electron",
            "table": self.table}
        with open(os.path.join(directory, "aspher.json"),'w') as out:
            out.write(json.dumps(model, indent=4))
        return True

    def run(self, directory, structure_file, hkl_file):
        if not self.save(directory):
            return
        subprocess.run(["./bin/discamb_ff.exe", structure_file, hkl_file], cwd=directory)

        
        