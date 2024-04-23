#!/usr/bin/env python
# -*- coding: UTF-8 -*-

# mod_lagna.py -- module lagna. All computations for lagna chart [D1 chart]
#
# Copyright (C) 2022 Shyam Bhat  <vicharavandana@gmail.com>
# Downloaded from "https://github.com/VicharaVandana/jyotishyam.git"
#
# This file is part of the "jyotishyam" Python library
# for computing Hindu jataka with sidereal lahiri ayanamsha technique
# using swiss ephemeries
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
Use Swiss ephemeris to calculate planetery position of 9 vedic astrology planets 
and lagna(Ascendant)
"""

from collections import namedtuple as struct
from operator import ge
import time
import swisseph as swe
import generic.mod_constants as c
# import generic.mod_general as self.gen

from mod_astrodata import AstroData
from generic.mod_general import General


class Lagna():

    data = AstroData()  # initate with birthdata
    gen = General()

    def __init__(self, name, gender, place, longitude, lattitude, timezone, year, month, day, hour, min, sec=0):
        self.data.birthdata = {"DOB": {"year": year,
                                       "month": month,
                                       "day": day
                                       },
                               "TOB": {"hour": hour,  # in 24 hour format
                                       "min": min,
                                       "sec": sec
                                       },
                               "POB": {"name": place,
                                       "lon": longitude,  # +ve for North and -ve for south
                                       "lat": lattitude,  # +ve for East and -ve for West
                                       "timezone": timezone
                                       },
                               "name": name,
                               "Gender": gender,
                               "Comments": "Yaaaa"
                               }
        self.bd = self.data.birthdata

    Date = struct("Date", ["year", "month", "day"])
    Place = struct("Place", ["latitude", "longitude", "timezone"])

    sidereal_year = 365.256360417   # From WolframAlpha

    # namah suryaya chandraya mangalaya ... rahuve ketuve namah
    swe.KETU = swe.PLUTO  # I've mapped Pluto to Ketu
    planet_list = [swe.SUN, swe.MOON, swe.MARS, swe.MERCURY, swe.JUPITER,
                   swe.VENUS, swe.SATURN, swe.MEAN_NODE,  # Rahu = MEAN_NODE
                   swe.KETU]

    def set_ayanamsa_mode(self): return swe.set_sid_mode(
        swe.SIDM_LAHIRI)  # Vedic astrology uses Lahiri ayanamsa
    # Default is FAGAN_BRADLEY ayanamsa in swiss ephemeries

    def reset_ayanamsa_mode(self): return swe.set_sid_mode(
        swe.SIDM_FAGAN_BRADLEY)

    ################################# FUNCTIONS #############################

    def get_planet_name(self, planet):
        names = {swe.SUN: "Sun", swe.MOON: "Moon", swe.MARS: "Mars",
                 swe.MERCURY: "Mercury", swe.JUPITER: "Jupiter", swe.VENUS: "Venus",
                 swe.SATURN: "Saturn", swe.MEAN_NODE: "Rahu", swe.KETU: "Ketu"}
        return names[planet]

    def get_planet_symbol(self, planet):
        symbols = {swe.SUN: "Su", swe.MOON: "Mo", swe.MARS: "Ma",
                   swe.MERCURY: "Me", swe.JUPITER: "Ju", swe.VENUS: "Ve",
                   swe.SATURN: "Sa", swe.MEAN_NODE: "Ra", swe.KETU: "Ke"}
        return symbols[planet]

    # Convert 23d 30' 30" to 23.508333 degrees

    def from_dms(self, degs, mins, secs): return degs + mins/60 + secs/3600

    # the inverse

    def to_dms_prec(self, deg):
        d = int(deg)
        mins = (deg - d) * 60
        m = int(mins)
        s = round((mins - m) * 60, 6)
        return [d, m, s]

    def to_dms(self, deg):
        d, m, s = self.to_dms_prec(deg)
        return [d, m, int(s)]

    # Make angle lie between [-180, 180) instead of [0, 360)

    def norm180(self, angle): return (angle - 360) if angle >= 180 else angle

    # Make angle lie between [0, 360)

    def norm360(self, angle): return angle % 360

    # Ketu is always 180° after Rahu, so same coordinates but different constellations
    # i.e if Rahu is in Pisces, Ketu is in Virgo etc

    def ketu(self, rahu): return (rahu + 180) % 360

    # Julian Day number as on (year, month, day) at 00:00 UTC

    def gregorian_to_jd(self, date): return swe.julday(
        date.year, date.month, date.day, 0.0)

    def jd_to_gregorian(self, jd): return swe.revjul(
        jd, swe.GREG_CAL)   # returns (y, m, d, h, min, s)

    get_nakshatra_name = ["Ashwini", "Bharani", "Kritika",
                          "Rohini", "Mrigashira", "Ardra",
                          "Punarvasu", "Pushya", "Ashlesha",
                          "Magha", "Purva Phalguni", "Uttara Phalguni",
                          "Hasta", "Chitra", "Swati",
                          "Vishaka", "Anurada", "Jyeshta",
                          "Mula", "Purva Ashadha", "Uttara Ashadha",
                          "Shravana", "Dhanishta", "Shatabhishak",
                          "Purva Bhadrapada", "Uttara Bhadrapada", "Revati"]

    def nakshatra_pada(self, longitude):
        """Gives nakshatra (0..26) and paada (1..4) in which given longitude lies"""
        # 27 nakshatras span 360°
        one_star = (360 / 27)  # = 13°20'
        # Each nakshatra has 4 padas, so 27 x 4 = 108 padas in 360°
        one_pada = (360 / 108)  # = 3°20'
        quotient = int(longitude / one_star)
        reminder = (longitude - quotient * one_star)
        pada = int(reminder / one_pada)
        # convert 0..26 to 1..27 and 0..3 to 1..4
        return [self.get_nakshatra_name[quotient], 1 + pada]

    def sidereal_longitude(self, jd, planet):
        """Computes nirayana (sidereal) longitude of given planet on jd"""
        self.set_ayanamsa_mode()
        (longi, myflags) = swe.calc_ut(jd, planet,
                                       flag=swe.FLG_SWIEPH | swe.FLG_SIDEREAL)
        self.reset_ayanamsa_mode()
        return self.norm360(longi[0])  # degrees

    def Is_Retrograde(self, jd, planet):
        """Checks if given planet is in retrograde motion on jd"""
        self.set_ayanamsa_mode()
        (longi, myflags) = swe.calc_ut(jd, planet,
                                       flag=swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_SIDEREAL)
        self.reset_ayanamsa_mode()
        return (longi[3] < 0)  # if speed is negative then its in retro

    def update_ascendant(self, jd, place):
        """Lagna (=ascendant) calculation at any given time & place
          It also updates most of lagna elements data, 
          except lagnesh_sign, lagnesh rashi and lagnesh dispositor"""
        lat, lon, tz = place
        jd_utc = jd - (tz / 24.)
        self.set_ayanamsa_mode()  # needed for swe.houses_ex()
        # returns two arrays, cusps and ascmc, where ascmc[0] = Ascendant
        nirayana_lagna = swe.houses_ex(
            jd_utc, lat, lon, flag=swe.FLG_SIDEREAL)[1][0]
        # 12 zodiac signs span 360°, so each one takes 30°
        # 0 = Mesha, 1 = Vrishabha, ..., 11 = Meena
        constellation = int(nirayana_lagna / 30)
        coordinates = self.to_dms(nirayana_lagna % 30)
        self.reset_ayanamsa_mode()
        # Updating the data from computed values
        # update position of ascendant
        self.data.lagna_ascendant["pos"]["deg"] = coordinates[0]
        self.data.lagna_ascendant["pos"]["min"] = coordinates[1]
        self.data.lagna_ascendant["pos"]["sec"] = coordinates[2]
        self.data.lagna_ascendant["pos"]["dec_deg"] = (nirayana_lagna % 30)

        # update nakshatra related self.data for ascendant
        nak_pad = self.nakshatra_pada(nirayana_lagna)
        self.data.lagna_ascendant["nakshatra"] = nak_pad[0]
        self.data.lagna_ascendant["pada"] = nak_pad[1]
        self.data.lagna_ascendant["nak-ruler"] = self.gen.ruler_of_nakshatra[nak_pad[0]]
        self.data.lagna_ascendant["nak-diety"] = self.gen.diety_of_nakshatra[nak_pad[0]]

        # update sign related self.data for ascendant
        self.data.lagna_ascendant["sign"] = self.gen.signs[constellation]
        self.data.lagna_ascendant["rashi"] = self.gen.rashis[constellation]
        self.data.lagna_ascendant["lagna-lord"] = self.gen.signlords[constellation]
        self.data.lagna_ascendant["sign-tatva"] = self.gen.signtatvas[constellation]

        # updating Status of Ascendant
        self.data.lagna_ascendant["status"] = c.PARTIAL

        return (1 + constellation)

    def update_planetaryData(self, jd, place):
        """Computes instantaneous planetary positions
          (i.e., which celestial object lies in which constellation)
          Also gives the nakshatra-pada division
          And updates the birth chart data for all the planets except those dependant on ascendant
        """
        jd_ut = jd - place.timezone / 24.

        for planet in self.planet_list:
            if planet != swe.KETU:
                nirayana_long = self.sidereal_longitude(jd_ut, planet)
                retro = self.Is_Retrograde(jd_ut, planet)
            else:  # Ketu
                # nirayana_long = ketu(self.sidereal_longitude(jd_ut, swe.RAHU))
                nirayana_long = self.ketu(
                    self.sidereal_longitude(jd_ut, swe.MEAN_NODE))
                retro = True  # ketu is always in retrograde

            # 12 zodiac signs span 360°, so each one takes 30°
            # 0 = Mesha, 1 = Vrishabha, ..., 11 = Meena
            constellation = int(nirayana_long / 30)
            coordinates = self.to_dms(nirayana_long % 30)

            # Update the data properly for the planet
            db_planet = self.data.lagna_planets[self.get_planet_name(
                planet)]  # get the proper planet container
            db_planet["retro"] = retro  # retrograde property

            # update position of the planet
            db_planet["pos"]["deg"] = coordinates[0]
            db_planet["pos"]["min"] = coordinates[1]
            db_planet["pos"]["sec"] = coordinates[2]
            db_planet["pos"]["dec_deg"] = (nirayana_long % 30)

            # update nakshatra related data for the planet
            nak_pad = self.nakshatra_pada(nirayana_long)
            db_planet["nakshatra"] = nak_pad[0]
            db_planet["pada"] = nak_pad[1]
            db_planet["nak-ruler"] = self.gen.ruler_of_nakshatra[nak_pad[0]]
            db_planet["nak-diety"] = self.gen.diety_of_nakshatra[nak_pad[0]]

            # update sign related data for the planet
            currentsign = self.gen.signs[constellation]
            db_planet["sign"] = currentsign
            db_planet["rashi"] = self.gen.rashis[constellation]
            dispositor = self.gen.signlords[constellation]
            db_planet["dispositor"] = dispositor
            db_planet["sign-tatva"] = self.gen.signtatvas[constellation]

            # Compute the house relation for the planet

            exhaltsign = self.gen.exhaltationSign_of_planet[db_planet["name"]]
            debilitsign = self.gen.debilitationSign_of_planet[db_planet["name"]]
            friends = db_planet["friends"]
            enemies = db_planet["enemies"]
            neutral = db_planet["nuetral"]

            if (currentsign == exhaltsign):  # first check for exhaltation
                db_planet["house-rel"] = c.EXHALTED
            elif (currentsign == debilitsign):  # next check for debilitated
                db_planet["house-rel"] = c.DEBILITATED
            elif (db_planet["name"] == dispositor):  # next check for own sign
                db_planet["house-rel"] = c.OWNSIGN
            elif (dispositor in friends):  # next check for friend sign
                db_planet["house-rel"] = c.FRIENDSIGN
            elif (dispositor in enemies):  # next check for enemy sign
                db_planet["house-rel"] = c.ENEMYSIGN
            elif (dispositor in neutral):  # next check for neutral sign
                db_planet["house-rel"] = c.NEUTRALSIGN
            else:
                db_planet["house-rel"] = "UNKNOWN"

            # updating Status of planet
            db_planet["status"] = c.PARTIAL

        return

    def compute_lagnaChart(self):
        birthday_julien = swe.julday(self.bd["DOB"]["year"],  # birth year
                                     self.bd["DOB"]["month"],  # birth month
                                     self.bd["DOB"]["day"],  # birth day
                                     # birth time in float
                                     ((self.bd["TOB"]["hour"]) + (self.bd["TOB"]["min"]) / \
                                      60. + (self.bd["TOB"]["sec"])/3600),
                                     )  # yyyy,mm,dd,time_24hr_format(hh + mm/60 + ss/3600)

        birth_place = self.Place(self.bd["POB"]["lon"],  # longitude
                                 self.bd["POB"]["lat"],  # lattitude
                                 self.bd["POB"]["timezone"]  # Timezone
                                 )
        # Compute ascendant related data
        lagna = self.update_ascendant(birthday_julien, birth_place)
        # Compute navagraha related data
        self.update_planetaryData(birthday_julien, birth_place)

        # update miscdata like maasa vaara tithi etc
        self.gen.update_miscdata(birthday_julien, birth_place,
                                 self.data.charts["user_details"])

        # computing benefics, malefics and neutral planets for given lagna
        self.gen.compute_BenMalNeu4lagna(
            lagna, self.data.D1["classifications"])

        # computing lagnesh related self.data for ascendant - not updated by update_ascendant()
        lagnesh = self.data.lagna_ascendant["lagna-lord"]  # get lagnesh
        # check the sign of lagnesh
        self.data.lagna_ascendant["lagnesh-sign"] = self.data.lagna_planets[lagnesh]["sign"]
        self.data.lagna_ascendant["lagnesh-rashi"] = self.data.lagna_planets[lagnesh]["rashi"]
        self.data.lagna_ascendant["lagnesh-disp"] = self.data.lagna_planets[lagnesh]["dispositor"]
        # updating Status of Ascendant
        self.data.lagna_ascendant["status"] = c.COMPUTED

        # computing house related self.data for planets - not updated by update_planetaryData()
        for planetname in self.data.lagna_planets:
            planet = self.data.lagna_planets[planetname]
            planet["house-num"] = self.gen.housediff(lagna,
                                                     self.gen.signnum(planet["sign"]))
            # updating Status of the planet
            planet["status"] = c.COMPUTED

        self.gen.update_houses(self.data.D1)

        # computing aspects and conjunction planets
        self.gen.compute_aspects(self.data.D1)
        self.gen.compute_aspectedby(self.data.D1)
        self.gen.compute_conjuncts(self.data.D1)

        # populating the classification part of divisional chart
        self.gen.populate_kendraplanets(self.data.D1)  # kendra planets
        self.gen.populate_trikonaplanets(self.data.D1)  # trikona planets
        self.gen.populate_trikplanets(self.data.D1)  # trik planets
        self.gen.populate_upachayaplanets(self.data.D1)  # upachaya planets
        self.gen.populate_dharmaplanets(self.data.D1)  # dharma planets
        self.gen.populate_arthaplanets(self.data.D1)  # artha planets
        self.gen.populate_kamaplanets(self.data.D1)  # kama planets
        self.gen.populate_mokshaplanets(self.data.D1)  # moksha planets

        return self.data.charts


# if __name__ == "__main__":
#     compute_lagnaChart()
