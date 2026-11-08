"""English-language GUI for Seestar photometry. Processing mathematics stays in seestar_photometry.py."""
from __future__ import annotations
import contextlib
import csv
import datetime
import configparser
import io
import json
import math
import os
import queue
import re
import sys
import threading
import traceback
from collections import Counter
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import pandas as pd
import seestar_photometry as engine

# Portable configuration: beside the launched EXE, or beside this .py file.
# Never use the temporary PyInstaller one-file extraction directory.
PROGRAM_DIR = (Path(sys.executable).resolve().parent if getattr(sys, 'frozen', False)
               else Path(__file__).resolve().parent)
SETTINGS = PROGRAM_DIR / 'settings.json'
LEGACY_APP_DIR = Path(os.environ.get('APPDATA', Path.home())) / 'SeestarPhotometry'
LEGACY_SETTINGS = LEGACY_APP_DIR / 'settings.json'
CARDS = LEGACY_APP_DIR / 'cards'  # Legacy cards remain readable; new cards go into the CSV folder


def load_site():
    # Prefer the portable settings. Migrate valid settings from the old AppData
    # location only when no local settings file already exists.
    for path in (SETTINGS, LEGACY_SETTINGS):
        try:
            d = validate_site(json.loads(path.read_text(encoding='utf-8')))
            if path == LEGACY_SETTINGS and not SETTINGS.exists():
                try:
                    SETTINGS.write_text(json.dumps(d, ensure_ascii=False, indent=2), encoding='utf-8')
                except OSError:
                    # Still allow reading the old settings; save will report any failure.
                    pass
            return d
        except (OSError, ValueError, KeyError, TypeError):
            continue
    return None


def validate_site(d):
    name = str(d['name']).strip()
    lat, lon, elev = (float(d[k]) for k in ('latitude','longitude','elevation'))
    if not name or not all(map(math.isfinite, (lat,lon,elev))) or not -90 <= lat <= 90 or not -180 <= lon <= 180:
        raise ValueError('Enter a name, latitude −90..90°, longitude −180..180° and finite elevation.')
    return {'name':name, 'latitude':lat, 'longitude':lon, 'elevation':elev}


def card_path(name, folder=None):
    if not name or re.search(r'[\\/:*?"<>|\x00-\x1f]', name) or name in ('.','..'):
        raise ValueError('Invalid variable star name for a Windows filename.')
    return (Path(folder) if folder is not None else CARDS) / (name + '.txt')


def read_card(path):
    cp = configparser.ConfigParser(interpolation=None)
    cp.read(path, encoding='utf-8-sig')
    if len(cp.sections()) != 1: raise ValueError('Card must contain exactly one [star] section.')
    return cp[cp.sections()[0]]


def write_card(fields, folder):
    name = fields['NAME'].strip()
    file = card_path(name, folder)
    cp = configparser.ConfigParser(interpolation=None)
    cp.optionxform = str
    cp[name] = {k:v for k,v in fields.items() if k != 'NAME' and str(v).strip()}
    Path(folder).mkdir(parents=True,exist_ok=True)
    with file.open('w', encoding='utf-8') as f: cp.write(f)
    return file


def fetch_vsx(ra, dec, radius=8.0):
    params=urlencode({'view':'api.list','ra':f'{ra:.8f}','dec':f'{dec:.8f}',
                      'radius':f'{radius/3600:.9f}','format':'json'})
    request=Request('https://www.aavso.org/vsx/index.php?'+params,
                    headers={'User-Agent':'SeestarPhotometry/2.0','Accept':'application/json'})
    with urlopen(request, timeout=15) as response:
        payload=json.loads(response.read().decode('utf-8-sig'))
    objects=payload.get('VSXObjects',payload)
    if isinstance(objects,dict): objects=objects.get('VSXObject',[])
    if isinstance(objects,dict): objects=[objects]
    found=[]
    for x in objects if isinstance(objects,list) else []:
        try:
            sep=float(engine.sep_arcsec(float(x['RA2000']),float(x['Declination2000']),ra,dec))
            if sep<=radius and str(x.get('Name','')).strip():
                found.append((sep,x))
        except (ValueError,TypeError,KeyError): pass
    return sorted(found,key=lambda e:e[0])


def get_vsx_field(obj, *names):
    """Read VSX fields independent of punctuation/casing in the API response."""
    def norm(key): return re.sub(r'[^a-z0-9]', '', str(key).lower())
    items = {norm(k): v for k, v in obj.items()}
    for key in names:
        val = items.get(norm(key))
        if isinstance(val, dict):
            val = val.get('#text', val.get('value', ''))
        if val is not None and str(val).strip() not in ('', '--', 'None'):
            return str(val).strip()
    return ''


def fetch_vsx_object(name):
    """Fetch a detailed VSX record, not just the abbreviated nearby-star listing."""
    params=urlencode({'view':'api.object','ident':name,'format':'json'})
    req=Request('https://www.aavso.org/vsx/index.php?'+params,
                headers={'User-Agent':'SeestarPhotometry/2.1','Accept':'application/json'})
    with urlopen(req, timeout=20) as response:
        payload=json.loads(response.read().decode('utf-8-sig'))
    objects=payload.get('VSXObjects',payload) if isinstance(payload,dict) else payload
    if isinstance(objects,dict): objects=objects.get('VSXObject',objects)
    if isinstance(objects,dict): objects=[objects]
    if not isinstance(objects,list) or not objects:
        raise ValueError('No detailed VSX record was returned for '+name)
    for obj in objects:
        if isinstance(obj,dict) and get_vsx_field(obj,'Name'):
            return obj
    raise ValueError('The detailed VSX response contained no star record.')


def vsx_card_fields(obj):
    """Only physical catalog parameters: no forecast minima from ephemeris tables."""
    aliases={
      'NAME':('Name',), 'RA':('RA2000','RA'), 'DEC':('Declination2000','Dec2000','Dec'),
      'TYPE':('VariabilityType','Type'), 'PERIOD':('Period',),
      'T0':('Epoch','EpochHJD','Epoch0','InitialEpoch'),
      'MAG_MIN':('MaxMag','MaximumMagnitude','Max','Maximum'),
      'MAG_MAX':('MinMag','MinimumMagnitude','Min','Minimum'),
      'ECLIPSE_DURATION':('RiseDuration','RiseEclipseDuration','EclipseDuration','Duration'),
    }
    fields={key:val for key,names in aliases.items() if (val:=get_vsx_field(obj,*names))}
    # VSX magnitudes can contain a passband suffix such as '9.71 V'.
    for key in ('MAG_MIN','MAG_MAX'):
        if key in fields:
            match=re.search(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)',fields[key])
            if match:fields[key]=match.group()
    if 'ECLIPSE_DURATION' in fields:
        match=re.search(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)',fields['ECLIPSE_DURATION'])
        if match:fields['ECLIPSE_DURATION']=match.group()
    return fields


def fetch_gaia(ra, dec, radius_arcsec=8.0):
    """Nearest Gaia DR3 J2016 position and proper motions via public Gaia TAP.

    Ambiguous astrometry is returned for review, not silently trusted.
    """
    radius=min(max(float(radius_arcsec), 0.25), 10.0)/3600.0
    query=("SELECT TOP 10 source_id,ra,dec,pmra,pmdec,phot_g_mean_mag "
           "FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS',ra,dec),"
           f"CIRCLE('ICRS',{ra:.9f},{dec:.9f},{radius:.10f})) "
           f"ORDER BY DISTANCE(POINT('ICRS',ra,dec),POINT('ICRS',{ra:.9f},{dec:.9f}))")
    params=urlencode({'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'csv','QUERY':query})
    request=Request('https://gea.esac.esa.int/tap-server/tap/sync?'+params,
                    headers={'User-Agent':'SeestarPhotometry/2.1'})
    with urlopen(request,timeout=30) as response:
        content=response.read().decode('utf-8-sig')
    rows=list(csv.DictReader(io.StringIO(content)))
    matches=[]
    for row in rows:
        try:
            if not row.get('pmra') or not row.get('pmdec'):continue
            sep=float(engine.sep_arcsec(float(row['ra']),float(row['dec']),ra,dec))
            matches.append((sep,row))
        except (ValueError,KeyError,TypeError):continue
    matches.sort(key=lambda item:item[0])
    return matches


class QueueWriter(io.TextIOBase):
    def __init__(self,messages,log_file):
        self.messages=messages
        self.log_file=log_file
    def write(self,s):
        if s:
            self.messages.put(('log',s))
            with self.log_file.open('a',encoding='utf-8') as f:f.write(s)
        return len(s)
    def flush(self): pass


class Application(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title('Seestar Photometry')
        self.geometry('900x740');self.minsize(740,560)
        self.messages=queue.Queue();self.running=False
        self.folder=tk.StringVar();self.star=tk.StringVar();self.radius=tk.StringVar(value='8.0')
        self.vsx_radius=tk.StringVar(value='8.0');self.min_points=tk.StringVar(value='100')
        self.show_median=tk.BooleanVar(value=True);self.mode=tk.StringVar(value='target')
        self.site=load_site()
        self.build_ui();self.change_mode()
        self.after(100,self.drain_messages)
        if self.site is None: self.after(350,lambda:self.edit_site(required=True))

    def build_ui(self):
        outer=ttk.Frame(self,padding=16);outer.pack(fill='both',expand=True)
        header=ttk.Frame(outer);header.pack(fill='x',pady=(0,12))
        ttk.Label(header,text='Seestar Photometry',font=('Segoe UI',18,'bold')).pack(side='left')
        ttk.Button(header,text='Observatory Settings',command=self.edit_site).pack(side='right')
        ttk.Button(header,text='Variable Cards',command=self.card_manager).pack(side='right',padx=8)
        self.site_label=ttk.Label(outer,text='',foreground='#555555');self.site_label.pack(anchor='w',pady=(0,10))
        self.update_site_label()
        ttk.Label(outer,text='Seestar CSV folder').pack(anchor='w')
        folder_row=ttk.Frame(outer);folder_row.pack(fill='x',pady=(5,12))
        ttk.Entry(folder_row,textvariable=self.folder).pack(side='left',fill='x',expand=True)
        ttk.Button(folder_row,text='Browse...',command=self.pick_folder).pack(side='left',padx=(8,0))
        ttk.Button(folder_row,text='Refresh',command=self.refresh_stars).pack(side='left',padx=(8,0))
        mode_row=ttk.LabelFrame(outer,text='Analysis mode',padding=10);mode_row.pack(fill='x',pady=(0,12))
        ttk.Radiobutton(mode_row,text='Target star — saved card / Time + Raw + Phase',variable=self.mode,
                        value='target',command=self.change_mode).pack(anchor='w')
        ttk.Radiobutton(mode_row,text='All Variables — Time plots only / VSX identification',variable=self.mode,
                        value='all',command=self.change_mode).pack(anchor='w')
        self.target_options=ttk.Frame(outer)
        ttk.Label(self.target_options,text='Target:').pack(side='left')
        self.star_select=ttk.Combobox(self.target_options,textvariable=self.star,state='readonly',width=30)
        self.star_select.pack(side='left',padx=(8,18))
        ttk.Label(self.target_options,text='Matching radius (″):').pack(side='left')
        ttk.Entry(self.target_options,textvariable=self.radius,width=8).pack(side='left',padx=8)
        self.vsx_options=ttk.Frame(outer)
        ttk.Label(self.vsx_options,text='VSX search radius (″):').pack(side='left')
        ttk.Entry(self.vsx_options,textvariable=self.vsx_radius,width=7).pack(side='left',padx=(8,15))
        ttk.Label(self.vsx_options,text='Minimum clean points:').pack(side='left')
        ttk.Entry(self.vsx_options,textvariable=self.min_points,width=7).pack(side='left',padx=(8,15))
        ttk.Checkbutton(self.vsx_options,text='Running median',variable=self.show_median).pack(side='left')
        self.note=ttk.Label(outer,text='Only CSV files in the selected folder are read (no subfolders). Results are saved in Results.',wraplength=750)
        self.note.pack(anchor='w',pady=(8,12))
        actions=ttk.Frame(outer);actions.pack(fill='x',pady=(0,12))
        self.start_btn=ttk.Button(actions,text='Run Analysis',command=self.start);self.start_btn.pack(side='left')
        ttk.Button(actions,text='Open Results Folder',command=self.open_results_folder).pack(side='left',padx=10)
        ttk.Label(outer,text='Processing log').pack(anchor='w',pady=(0,5))
        logframe=ttk.Frame(outer);logframe.pack(fill='both',expand=True)
        self.log=tk.Text(logframe,wrap='word',state='disabled',height=12,font=('Consolas',9))
        self.log.pack(side='left',fill='both',expand=True)
        bar=ttk.Scrollbar(logframe,command=self.log.yview);bar.pack(side='right',fill='y')
        self.log.configure(yscrollcommand=bar.set)
        ttk.Separator(outer).pack(fill='x',pady=(12,6))
        ttk.Label(outer,text='Zsolt Szamosvári  |  zsolt.szamosv@gmail.com',anchor='center',foreground='#666666').pack(fill='x')

    def update_site_label(self):
        self.site_label.configure(text=('Observatory: '+self.site['name']+f"  ({self.site['latitude']:.6f}°, {self.site['longitude']:.6f}°)"
                                        if self.site else 'Observatory: NOT CONFIGURED — set your location before processing'))

    def edit_site(self,required=False):
        win=tk.Toplevel(self);win.title('Observatory Settings');win.geometry('465x320');win.transient(self);win.grab_set()
        d=self.site or {}
        names=[('name','Observatory name'),('latitude','Latitude (degrees; north +)'),
               ('longitude','Longitude (degrees; east +)'),('elevation','Elevation above sea level (m)')]
        fields={}
        frame=ttk.Frame(win,padding=18);frame.pack(fill='both',expand=True)
        ttk.Label(frame,text='HJD correction uses these coordinates for every new run.',wraplength=420).pack(anchor='w',pady=(0,12))
        for key,label in names:
            row=ttk.Frame(frame);row.pack(fill='x',pady=4)
            ttk.Label(row,text=label,width=30).pack(side='left')
            var=tk.StringVar(value=str(d.get(key,'')));ttk.Entry(row,textvariable=var,width=20).pack(side='left',fill='x',expand=True)
            fields[key]=var
        ttk.Label(frame,text='Important: use the observing location, not your time zone.',foreground='#555555').pack(anchor='w',pady=12)
        def save():
            try:
                site=validate_site({k:v.get().replace(',','.') if k!='name' else v.get() for k,v in fields.items()})
                SETTINGS.write_text(json.dumps(site,ensure_ascii=False,indent=2),encoding='utf-8')
                self.site=site;self.update_site_label();win.destroy()
            except (OSError,ValueError,KeyError) as e:messagebox.showerror('Invalid observatory settings',str(e),parent=win)
        ttk.Button(frame,text='Save Observatory',command=save).pack(anchor='e')
        if required: win.protocol('WM_DELETE_WINDOW',lambda:win.destroy())

    def pick_folder(self):
        p=filedialog.askdirectory(title='Choose Seestar CSV folder')
        if p:self.folder.set(p);self.refresh_stars()

    def refresh_stars(self):
        folder=Path(self.folder.get().strip())
        if not folder.is_dir():messagebox.showwarning('Folder','Select an existing folder.');return
        stars={}
        for base in (CARDS,folder):
            if not base.is_dir():continue
            for path in base.glob('*.txt'):
                if path.name.lower()=='readme.txt':continue
                try:
                    cp=engine.read_targets(path)
                    if len(cp.sections())==1 and path.stem.casefold()==cp.sections()[0].casefold():
                        stars[cp.sections()[0].casefold()]=cp.sections()[0]
                except Exception:pass
        names=sorted(stars.values(),key=str.casefold)
        old=self.star.get()
        self.star_select['values']=names
        self.star.set(old if old in names else names[0] if names else '')
        if not names and self.mode.get()=='target':
            messagebox.showinfo('No variable cards','Create or import a variable card using Variable Cards.')

    def change_mode(self):
        self.target_options.pack_forget();self.vsx_options.pack_forget()
        if self.mode.get()=='all':self.vsx_options.pack(fill='x',pady=(0,12),before=self.note)
        else:self.target_options.pack(fill='x',pady=(0,12),before=self.note)

    def card_manager(self):
        win=tk.Toplevel(self);win.title('Variable Card Manager');win.geometry('740x730');win.minsize(650,590);win.transient(self)
        frame=ttk.Frame(win,padding=12);frame.pack(fill='both',expand=True)
        top=ttk.Frame(frame);top.pack(fill='x',pady=(0,10))
        ttk.Label(top,text='Saved / existing cards:').pack(side='left')
        choice=tk.StringVar(); combo=ttk.Combobox(top,textvariable=choice,state='readonly',width=30);combo.pack(side='left',padx=8)
        props=[('NAME','Variable ID'),('GAIA_DR3','Gaia DR3 ID (optional)'),('RA','RA (degrees)'),('DEC','Dec (degrees)'),
               ('COORD_EPOCH','Coordinate epoch (year)'),('PMRA','PM RA (mas/year)'),('PMDEC','PM Dec (mas/year)'),
               ('T0','T0 — HJD (verify manually)'),('PERIOD','Period (days)'),('TYPE','Variability type'),
               ('MAG_MIN','Bright magnitude'),('MAG_MAX','Faint magnitude'),('ECLIPSE_DURATION','Eclipse duration (% of P)'),
               ('SOURCE','Source / reference')]
        values={}; panel=ttk.Frame(frame);panel.pack(fill='x')
        for k,label in props:
            line=ttk.Frame(panel);line.pack(fill='x',pady=2)
            ttk.Label(line,text=label,width=32).pack(side='left')
            var=tk.StringVar(value='');ent=ttk.Entry(line,textvariable=var,width=42);ent.pack(side='left',fill='x',expand=True)
            values[k]=var
        verified=tk.BooleanVar(value=False)
        ttk.Checkbutton(frame,text='I checked T0, its time scale, period and reference epoch.',variable=verified).pack(anchor='w',pady=(10,4))
        ttk.Label(frame,text='Enter a VSX name and use VSX Search to fill catalog data, or Suggest from CSV first.\nT0 from VSX is a PROPOSAL only: confirm manually before saving. Forecast minima are not used.',
                  wraplength=680,foreground='#555555').pack(anchor='w',pady=(0,8))
        def set_card(v):
            for k,var in values.items():var.set(str(v.get(k,'')))
            verified.set(str(v.get('T0_VERIFIED','')).upper()=='YES')
        def cards_list():
            items={}
            for base in (CARDS,Path(self.folder.get().strip()) if self.folder.get().strip() else Path('/nonexistent')):
                if base.is_dir():
                    for path in base.glob('*.txt'):
                        try:
                            card=read_card(path)
                            if path.stem.casefold()==card.name.casefold():items[card.name]=path
                        except Exception:pass
            return items
        paths={}
        def reload():
            nonlocal paths
            paths=cards_list();combo['values']=sorted(paths,key=str.casefold)
        def select(event=None):
            name=choice.get()
            if name in paths:
                card=read_card(paths[name])
                # ConfigParser stores all keys lowercase in dict(); normalize.
                s={k.upper():v for k,v in card.items()};s['NAME']=card.name;set_card(s)
        combo.bind('<<ComboboxSelected>>',select)
        def new():
            choice.set('');set_card({})
            for k,v in [('SOURCE','VSX / user verified')]:values[k].set(v)
        def save():
            data={k:v.get().strip() for k,v in values.items()}
            try:
                folder=Path(self.folder.get().strip())
                if not folder.is_dir():raise ValueError('Choose the CSV source folder in the main window first.')
                if not verified.get():raise ValueError('Please verify the ephemeris values first.')
                for k in ('RA','DEC','COORD_EPOCH','PMRA','PMDEC','T0','PERIOD'):
                    if not data[k]:raise ValueError(f'{k} is required.')
                    n=float(data[k]);
                    if not math.isfinite(n):raise ValueError(f'{k} must be finite.')
                if float(data['PERIOD'])<=0 or not 0<=float(data['RA'])<360 or not -90<=float(data['DEC'])<=90:
                    raise ValueError('Invalid period or coordinates.')
                data['T0_TYPE']='HJD';data['T0_VERIFIED']='YES'
                dest=card_path(data['NAME'],folder)
                if dest.exists() and not messagebox.askyesno('Overwrite card',f'Overwrite {dest.name}?',parent=win):return
                write_card(data,folder);reload();choice.set(data['NAME']);self.refresh_stars()
                messagebox.showinfo('Saved',f'Card saved: {dest}',parent=win)
            except (ValueError,OSError) as e:messagebox.showerror('Invalid card',str(e),parent=win)
        def vsx_fetch():
            try:
                name=values['NAME'].get().strip()
                ra_text=values['RA'].get().strip();dec_text=values['DEC'].get().strip()
                ra=dec=None
                if ra_text and dec_text:
                    ra=float(ra_text.replace(',','.'));dec=float(dec_text.replace(',','.'))
                radius=float(self.vsx_radius.get().replace(',','.'))
                if name:
                    detailed=fetch_vsx_object(name)
                else:
                    if ra is None or dec is None:
                        raise ValueError('Enter the variable name, or first use Suggest from CSV.')
                    matches=fetch_vsx(ra,dec,radius)
                    if not matches:raise ValueError('No VSX identification near the supplied coordinates.')
                    if len(matches)>1:
                        raise ValueError('Ambiguous VSX matches; enter a unique variable name:\n'+
                            '\n'.join(f"{x.get('Name')} — {sep:.2f}″" for sep,x in matches[:8]))
                    detailed=fetch_vsx_object(matches[0][1]['Name'])
                extracted=vsx_card_fields(detailed)
                if not extracted.get('NAME'):raise ValueError('VSX did not provide a variable name.')
                # Keep original Gaia astrometry rather than overwrite with VSX J2000 positions.
                star_ra=extracted.get('RA');star_dec=extracted.get('DEC')
                gaia_note=''
                if star_ra and star_dec:
                    try:
                        matches=fetch_gaia(float(star_ra),float(star_dec),radius)
                        if len(matches)==1 or (len(matches)>1 and matches[1][0]>max(2.0,2.5*matches[0][0])):
                            sep,row=matches[0]
                            extracted.update(GAIA_DR3=row['source_id'], RA=row['ra'], DEC=row['dec'],
                                             COORD_EPOCH='2016.0', PMRA=row['pmra'], PMDEC=row['pmdec'])
                            gaia_note=f'Gaia DR3 astrometry loaded ({sep:.2f}″).'
                        elif matches:
                            gaia_note='WARNING: multiple nearby Gaia matches; astrometry needs manual identification.'
                            for key in ('RA','DEC'):extracted.pop(key,None)
                        else:
                            gaia_note='No Gaia DR3 source with proper motion found. Complete Gaia details before saving.'
                            for key in ('RA','DEC'):extracted.pop(key,None)
                    except Exception as exc:
                        gaia_note=f'Gaia service unavailable ({type(exc).__name__}: {exc}). Astrometry needs review.'
                        for key in ('RA','DEC'):extracted.pop(key,None)
                for key,val in extracted.items():
                    if key in values:values[key].set(str(val))
                values['SOURCE'].set('AAVSO VSX + Gaia DR3 (verify T0)')
                verified.set(False)
                missing=[k for k in ('GAIA_DR3','RA','DEC','COORD_EPOCH','PMRA','PMDEC','T0','PERIOD') if not values[k].get().strip()]
                messagebox.showinfo('Catalog data loaded',
                    f"VSX: {extracted['NAME']}\n{gaia_note}\n"
                    + ('Missing fields: '+', '.join(missing)+'\n' if missing else '')
                    + 'The T0 value has NOT been approved. Confirm it before saving.',parent=win)
            except Exception as e:
                messagebox.showerror('VSX / Gaia retrieval failed',str(e),parent=win)
        # Keep candidates in their CSV row order; each row is a different star.
        # Do NOT collapse them into the most frequently repeated target position.
        csv_candidates=[]
        candidate_index=-1
        candidate_note=tk.StringVar(value='No CSV candidate selected')

        def load_csv_candidates():
            folder=Path(self.folder.get().strip())
            if not folder.is_dir():
                raise ValueError('Select a CSV source folder in the main window first.')
            files=sorted(folder.glob('*.csv'))
            if not files:
                raise ValueError('No CSV files found in the selected folder.')
            # A single exposure already contains the complete list of stars.
            # Using its row order avoids reading thousands of night-long CSVs.
            for file in files:
                try:
                    d=pd.read_csv(file,usecols=lambda c:c in ('Name','RA','DEC','Is Variable','Mag(obs)'),low_memory=False)
                    if not {'RA','DEC','Is Variable'}.issubset(d.columns):continue
                    isvar=d['Is Variable'].astype(str).str.strip().str.lower().isin(('true','1','yes'))
                    candidates=[]
                    for rowno,row in d.loc[isvar].iterrows():
                        ra=pd.to_numeric(row['RA'],errors='coerce')
                        dec=pd.to_numeric(row['DEC'],errors='coerce')
                        if not math.isfinite(ra) or not math.isfinite(dec) or not (0<=ra<360 and -90<=dec<=90):continue
                        name=str(row.get('Name','')).strip()
                        if name.lower() in ('nan','none','null'):name=''
                        mag=pd.to_numeric(row.get('Mag(obs)',float('nan')),errors='coerce')
                        candidates.append({'ra':float(ra),'dec':float(dec),'name':name,
                                           'mag':float(mag),'row':int(rowno)+2})
                    if candidates:
                        return candidates,file.name
                except (ValueError,KeyError,OSError):
                    continue
            raise ValueError('No valid Is Variable=True rows found in the CSV files.')

        def show_candidate(index):
            nonlocal candidate_index
            if not csv_candidates:return
            candidate_index=index%len(csv_candidates)
            obj=csv_candidates[candidate_index]
            new()
            values['RA'].set(f"{obj['ra']:.9f}")
            values['DEC'].set(f"{obj['dec']:.9f}")
            if obj['name']:values['NAME'].set(obj['name'])
            verified.set(False)
            candidate_note.set(f"Candidate {candidate_index+1}/{len(csv_candidates)}: "
                               f"{obj['name'] or 'Unnamed'} | RA {obj['ra']:.6f}°, Dec {obj['dec']:+.6f}°")

        def from_csv():
            nonlocal csv_candidates,candidate_index
            try:
                csv_candidates,filename=load_csv_candidates()
                show_candidate(0)
                messagebox.showinfo('CSV candidates',
                    f'{len(csv_candidates)} stars marked Is Variable=True found in:\n{filename}\n\n'
                    'Use Next Candidate to visit each one in its CSV order, then VSX Search to identify it.\n'
                    'Only manually verified cards can be saved.',parent=win)
            except Exception as exc:
                messagebox.showerror('CSV candidates',str(exc),parent=win)

        def next_candidate():
            if not csv_candidates:
                from_csv()
            else:
                show_candidate(candidate_index+1)
        buttons=ttk.Frame(frame);buttons.pack(fill='x',pady=6)
        ttk.Button(buttons,text='New Card',command=new).pack(side='left')
        ttk.Button(buttons,text='Suggest from CSV',command=from_csv).pack(side='left',padx=5)
        ttk.Button(buttons,text='Next Candidate',command=next_candidate).pack(side='left',padx=5)
        ttk.Button(buttons,text='VSX Search',command=vsx_fetch).pack(side='left',padx=5)
        ttk.Button(buttons,text='Save Verified Card',command=save).pack(side='right')
        ttk.Label(frame,textvariable=candidate_note,wraplength=680,foreground='#555555').pack(anchor='w',pady=(2,4))
        reload();new()

    def open_results_folder(self):
        folder=Path(self.folder.get().strip())
        output=folder/'Results'
        if not output.is_dir():messagebox.showinfo('No results','The Results directory does not exist yet.');return
        try:os.startfile(str(output))
        except OSError as e:messagebox.showerror('Cannot open folder',str(e))

    def append(self,s):
        self.log.configure(state='normal');self.log.insert('end',s);self.log.see('end');self.log.configure(state='disabled')

    def start(self):
        folder=Path(self.folder.get().strip())
        if not folder.is_dir():messagebox.showerror('Folder','Select an existing CSV folder.');return
        if not self.site:messagebox.showerror('Observatory','Set your observatory location first.');self.edit_site(required=True);return
        mode=self.mode.get();star=self.star.get().strip()
        if mode=='target' and not star:messagebox.showerror('Target','Select or create a variable card.');return
        try:
            radius=float(self.radius.get().replace(',','.'))
            vsx_radius=float(self.vsx_radius.get().replace(',','.'))
            points=int(self.min_points.get())
            if not 0<radius<=3600 or not 0<vsx_radius<=60 or points<2:raise ValueError
        except ValueError:messagebox.showerror('Parameters','Check matching radius, VSX radius and minimum points.');return
        # Prepare selected card in input directory for unchanged engine main().
        if mode=='target':
            local=folder/(star+'.txt');library=card_path(star)
            if local.exists():source=local
            elif library.is_file():source=library
            else:messagebox.showerror('Missing card','Target card not found.');return
            try:
                cp=engine.read_targets(source);t=engine.get_target(cp,star)
                if cp[star].get('T0_VERIFIED','').upper() not in ('YES',''):
                    messagebox.showwarning('Unverified card','Please verify T0 in Variable Cards first.');return
                # Legacy cards retain their existing status; never alter them silently.
            except Exception as e:messagebox.showerror('Invalid card',str(e));return
        else:source=None
        self.running=True;self.start_btn.configure(state='disabled')
        self.append(f'\n=== {mode} | {star if mode=="target" else "All Variables"} | {folder} ===\n')
        threading.Thread(target=self.process,args=(folder,star,radius,mode,vsx_radius,points,self.show_median.get(),source,dict(self.site)),daemon=True).start()

    def process(self,folder,star,radius,mode,vsx_radius,points,median,source,site):
        log_file=folder/'Results'/'photometry_processing.log'
        log_file.parent.mkdir(parents=True,exist_ok=True)
        writer=QueueWriter(self.messages,log_file)
        writer.write('\n=== '+datetime.datetime.now().isoformat(timespec='seconds')+
                     f' | {mode} | {star} | {folder} ===\n')
        try:
            with contextlib.redirect_stdout(writer),contextlib.redirect_stderr(writer):
                if mode=='all':
                    engine.all_variables(folder,match_arcsec=2.0,vsx_radius_arcsec=vsx_radius,
                                         min_points=points,show_median=median,site_name=site['name'])
                else:
                    # Temporarily stage the selected card for this run only if needed.
                    dest=folder/(star+'.txt');existing=dest.read_bytes() if dest.exists() else None
                    try:
                        if source.resolve()!=dest.resolve():dest.write_bytes(source.read_bytes())
                        # Avoid modifying global argv in the GUI worker thread: pass a local
                        # argument vector only during the synchronous engine.main() call.
                        previous=sys.argv
                        try:
                            sys.argv=['seestar_photometry.py','--folder',str(folder),'--star',star,
                                      '--max-sep',str(radius),'--site-settings',str(SETTINGS)]
                            engine.main()
                        finally:sys.argv=previous
                    finally:
                        if source.resolve()!=dest.resolve():
                            if existing is None:dest.unlink(missing_ok=True)
                            else:dest.write_bytes(existing)
            self.messages.put(('done',True))
        except BaseException:
            writer.write('\nERROR:\n'+traceback.format_exc()+'\n')
            self.messages.put(('done',False))

    def drain_messages(self):
        try:
            while True:
                action,value=self.messages.get_nowait()
                if action=='log':self.append(value)
                elif action=='done':
                    self.running=False;self.start_btn.configure(state='normal')
                    if value:self.append('\nProcessing complete.\n');messagebox.showinfo('Complete','Results saved to the Results folder.')
                    else:messagebox.showerror('Processing error', 'Processing failed. Log: '+str(Path(self.folder.get().strip())/'Results'/'photometry_processing.log'))
        except queue.Empty:pass
        self.after(100,self.drain_messages)


if __name__=='__main__':Application().mainloop()
