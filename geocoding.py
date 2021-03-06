#!/usr/bin/env python
# coding: utf-8

import os
import pandas as pd
import numpy as np 
from mapbox import Geocoder
import json

dataset = pd.read_csv("elenco_esercizi_commercio_in_sede_fissa_anno_2018.csv", sep = ';', header='infer', encoding='latin-1')

dataset['quartiere_settore'] = dataset.ESERCIZIO_VIA+'  '+dataset.ESERCIZIO_CIVICO+' '+dataset.QUARTIERE+' Bologna'
dataset['lat'] = float
dataset['lon'] = float

token = yourtoken
geocoder = Geocoder(access_token=token)

for i in range(len(dataset)):
    try:
        response = geocoder.forward(dataset.quartiere_settore[i])
        response = response.content
        output = json.loads(response)
        coordinates = output["features"][0]['geometry']['coordinates']
        dataset.iat[i,29] = coordinates[1]
        dataset.iat[i,30] = coordinates[0]
    except:
        dataset.iat[i,29] = np.nan
        dataset.iat[i,30] = np.nan

dataset.to_excel('geocoded.xlsx')

