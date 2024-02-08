#!/usr/bin/env python
# coding: utf-8

import pygrib, pytz, subprocess, sys
#import json, math, matplotlib, os, pygrib, pytz, struct, subprocess, sys
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import metpy.calc as mpcalc
import numpy as np
import pandas as pd
import pandas.tseries.offsets as offsets
from datetime import datetime, timedelta
#from itertools import repeat
from metpy.units import units
#from scipy.interpolate import griddata
#from scipy.ndimage import gaussian_filter, maximum_filter, minimum_filter
#from urllib.request import urlopen
#from scipy import interpolate
#from scipy.interpolate import interp2d, RectBivariateSpline, RegularGridInterpolator
#import argparse

markersize_0 = 2 # マーカーサイズ
char_size=16 # 文字サイズ
barb_length=8 # 矢羽の長さ
dlon,dlat=5,5   # 緯線・経線は1度ごと

# 描画地域と描画時刻の設定
if len(sys.argv) == 2:
    year = arg[:4]
    month = arg[4:6]
    day = arg[6:8]
    hour = arg[8:10]
else:
    jst = pytz.timezone('Asia/Tokyo')
    dt = datetime.now(jst) - timedelta(minutes=30)
    year=dt.year
    month=dt.month
    day=dt.day
    hour=dt.hour

dt = datetime(int(year), int(month), int(day), int(hour), 0)

# 初期値から5時間後に配信される
modeltime = pd.Timestamp(year,month,day,hour,0) - offsets.Hour(9) - offsets.Hour(5)
# MSMは03シリーズ
base_time = modeltime.replace(hour=modeltime.hour - (modeltime.hour % 3), minute=0, second=0)  

# データファイル名の指定
day_dir = base_time.strftime("%Y/%m/%d")
basename = "Z__C_RJTD_{}00_MSM_GPV_Rjp_L-pall_FH00-15_grib2.bin".format(base_time.strftime("%Y%m%d%H%M"))
url      = "http://database.rish.kyoto-u.ac.jp/arch/jmadata/data/gpv/original/{}/{}".format(day_dir, basename)
subprocess.run("wget {} -P ./".format(url), shell=True)       
grbs = pygrib.open(basename)

# 配列の宣言
grid_lon, grid_lat = np.meshgrid(np.arange(120, 150 + 0.0625, 0.125), np.arange(22.4, 47.6, 0.1))

# データを取得する
height = np.flip(grbs.select(parameterName='Geopotential height', level=500, forecastTime=3)[0].data()[0], axis=0)
tmp = np.flip(grbs.select(parameterName='Temperature', level=500, forecastTime=3)[0].data()[0] - 273.15, axis=0)
u = np.flip(grbs.select(parameterName='u-component of wind', level=850, forecastTime=3)[0].data()[0], axis=0)
v = np.flip(grbs.select(parameterName='v-component of wind', level=850, forecastTime=3)[0].data()[0], axis=0)
# 相当温位
tmp850 = np.flip(grbs.select(parameterName='Temperature', level=850, forecastTime=3)[0].data()[0] * units('K'), axis=0)
rh = np.flip(grbs.select(parameterName='Relative humidity', level=850, forecastTime=3)[0].data()[0] / 100, axis=0)
dewpoint = mpcalc.dewpoint_from_relative_humidity(tmp850, rh)
ept = mpcalc.equivalent_potential_temperature(850*units('hPa'), tmp850, dewpoint)

height = gaussian_filter(height, sigma=2.0)
tmp = gaussian_filter(tmp, sigma=2.0)
ept = gaussian_filter(ept, sigma=2.0)

# 図法指定                                                                             
proj = ccrs.PlateCarree()
# 図のSIZE指定inch                                                                        
fig = plt.figure(figsize=(8,6))
# 余白設定                                                                                
plt.subplots_adjust(left=0.04, right=1.1, bottom=0.0, top=1.0)                  
# 作図                                                                                    
ax = fig.add_subplot(1, 1, 1, projection=proj)
ax.set_extent([120, 150, 22.4, 47.6], proj)

# グリッド線を引く                                                               
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, alpha=0.8)
gl.xlocator = mticker.FixedLocator(np.arange(-180,180,5))
gl.ylocator = mticker.FixedLocator(np.arange(-90,90,5))



# プロット
cont = plt.contour(grid_lon, grid_lat, tmp, levels=np.arange(-60, 45, 3), linewidths=2, linestyles='solid', colors='red')
plt.clabel(cont, fontsize=30)

cont = plt.contour(grid_lon, grid_lat, height, levels=np.arange(5100, 6000, 60), linewidths=2, colors='black')
plt.clabel(cont, fontsize=30)

# 等圧線をプロット
#cont = plt.contour(grid_lon, grid_lat, ept, levels=np.arange(240, 390, 3), linewidths=2, colors='green')
plt.contourf(grid_lon, grid_lat, ept, cmap='coolwarm', levels=np.arange(246, 381, 9))
cb = plt.colorbar(orientation="vertical", shrink=0.6)    
#plt.clabel(cont, fontsize=30)

# ベクトルの間引き間隔
stride = 5

# データを間引く
grid_lon_sparse = grid_lon[::stride, ::stride]
grid_lat_sparse = grid_lat[::stride, ::stride]
u_sparse = u[::stride, ::stride]
v_sparse = v[::stride, ::stride]

ax.barbs(grid_lon_sparse, grid_lat_sparse, u_sparse, v_sparse, length=4, transform=latlon_proj)

# 等圧線のラベルを付ける
plt.clabel(cont, fontsize=30)

# 海岸線
ax.coastlines(resolution='10m', linewidth=1.6, color='black')  
            
# 図の説明
plt.title('{}'.format("AMeDAS and RadarGPV"), loc='left',size=20)
plt.title('Valid Time: {}'.format(dt), loc='right',size=20);
#plt.savefig("{}.jpg".format(time.strftime("%Y%m%d%H%M")), format="jpg")
plt.savefig("latest_upper.jpg", format="jpg")
