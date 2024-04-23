from chart_calc.mod_lagna import Lagna


# data = AstroData()
mod_lagna = Lagna(
    name='Kuldeep',
    day=18,
    month=7,
    year=2001,
    gender='Male',
    hour=5,
    min=30,
    lattitude=70.3,
    longitude=34.3,
    place='Bhiwani',
    timezone=+5.3,
)
da = mod_lagna.compute_lagnaChart()
print(da)
