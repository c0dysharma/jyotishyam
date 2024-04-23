import generic.mod_constants as c
from mod_astrodata import AstroData
import drawCharts.mod_drawChart as dc
import json


class JsonOut:
    data = AstroData()

    def dump_astrodata_injson(self, ):
        with open('./database/astrodata.json', 'w') as json_astrodatafile:
            json.dump(dict(self.data.charts), json_astrodatafile, indent=4)
        return

    def load_birthdatas(self, ):
        with open('./database/birthdatas.json', 'r') as json_birthfile:
            self.data.birthdatas = json.loads(json_birthfile.read())
        return

    def load_drawChartConfig(self, ):
        with open('./drawCharts/chartDraw_cfg.json', 'r') as json_birthfile:
            dc.chartCfg = json.loads(json_birthfile.read())
        return

    def get_birthdata(self, id):
        self.load_birthdatas()
        needed_birthdata = self.data.birthdatas.get(id, "NOT_FOUND")
        # check if the element with ID exists
        if ("NOT_FOUND" == needed_birthdata):
            # Element cannot be fetched
            print(f'Given ID {id} doesnt exist. So can not be fetched')
        return needed_birthdata

    def dump_birthdatas_injson(self, ):
        with open('./database/birthdatas.json', 'w') as json_birthfile:
            json.dump(dict(self.data.birthdatas), json_birthfile, indent=4)
        return

    def add_birthdata2DB(self, birthdata, id):
        # check if the element with ID already exists
        if ("NOT_FOUND" == self.data.birthdatas.get(id, "NOT_FOUND")):
            # New element -can be added
            self.data.birthdatas[id] = birthdata
            return True
        else:  # element already exist and so not possible to add
            print(f'Given ID {id} already exists. So can not be added')
            return False
