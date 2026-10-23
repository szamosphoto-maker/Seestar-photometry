"""Seestar CSV photometry -> photometry-quality report, time and phase plots.

Usage (Windows):
    python seestar_photometry.py
    python seestar_photometry.py --folder "D:/Astro/MyCSV"

Place one or more star-name TXT files (e.g. VW Cep.txt) in the CSV folder. Install:
    python -m pip install numpy pandas matplotlib astropy

HJD phase calculation requires Astropy; without it the time plots and tables are
still generated, and the phase plot is skipped, rather than mixing JD and HJD.
"""
from __future__ import annotations
import argparse
import configparser
import io
import zipfile
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

REQ = {'RA', 'DEC', 'Mag(obs)', 'Mag(err)', 'EXP start(JD)'}

def sep_arcsec(ra, dec, ra0, dec0):
    dra = ((ra-ra0+180) % 360)-180
    return 3600*np.hypot(dra*np.cos(np.deg2rad(dec0)),dec-dec0)

def read_targets(path):
    cp=configparser.ConfigParser(interpolation=None)
    if not cp.read(path,encoding='utf-8-sig'):
        raise FileNotFoundError(f'Nem található a katalógus: {path}')
    return cp

def get_target(cp,name):
    if name not in cp: raise ValueError(f'Ismeretlen cél: {name}. Választható: {", ".join(cp.sections())}')
    s=cp[name]
    target = {k:float(s[k]) for k in ('RA','DEC','COORD_EPOCH','PMRA','PMDEC','T0','PERIOD')}
    target.update({'T0_TYPE':s.get('T0_TYPE','HJD').upper(), 'NAME':name,
                   'GAIA_DR3':s.get('GAIA_DR3',''), 'SITE':s.get('SITE','').strip()})
    # Coordinates are optional, but when provided they must be complete.
    keys = ('LATITUDE','LONGITUDE','ELEVATION')
    provided = [bool(s.get(k,'').strip()) for k in keys]
    if any(provided) and not all(provided):
        raise ValueError('Az észlelőhelyhez LATITUDE, LONGITUDE és ELEVATION is szükséges.')
    if all(provided):
        for k in keys:
            target[k] = float(s[k])
            if not np.isfinite(target[k]):
                raise ValueError(f'Érvénytelen észlelőhely-adat: {k}')
        if not -90 <= target['LATITUDE'] <= 90:
            raise ValueError('LATITUDE értéke -90 és +90 fok közötti legyen.')
        if not -180 <= target['LONGITUDE'] <= 180:
            raise ValueError('LONGITUDE értéke -180 és +180 fok közötti legyen.')
        target['LOCATION_MODE'] = 'observatory'
    else:
        target['LOCATION_MODE'] = 'geocenter'
    return target

def predicted_position(t, year):
    dt=year-t['COORD_EPOCH']
    return (t['RA']+t['PMRA']*dt/(3.6e6*np.cos(np.deg2rad(t['DEC']))),
            t['DEC']+t['PMDEC']*dt/3.6e6)

def collect(files, t, threshold=8.0):
    """Read single-exposure and multi-exposure Seestar CSVs into one time series.

    A single-exposure file has (effectively) one exposure JD shared by its
    numerous stars. A tracked-star export has a distinct JD per exposure.
    """
    rows, errors = [], []
    seen = unmatched = 0
    multi_files = multi_rows = 0
    for filename, raw in files:
        seen += 1
        try:
            df = pd.read_csv(io.BytesIO(raw), low_memory=False)
            if not REQ.issubset(df.columns):
                errors.append(f'{filename}: hiányzó kötelező oszlop')
                continue
            jds = pd.to_numeric(df['EXP start(JD)'], errors='coerce').to_numpy(float)
            finite_jd = np.isfinite(jds)
            if not finite_jd.any():
                errors.append(f'{filename}: nincs érvényes JD')
                continue
            # 1 second tolerance avoids numerical JD rounding causing false detection.
            is_multi = (np.nanmax(jds) - np.nanmin(jds)) * 86400 > 1.0
            if is_multi:
                multi_files += 1
            ra = pd.to_numeric(df.RA, errors='coerce').to_numpy(float)
            dec = pd.to_numeric(df.DEC, errors='coerce').to_numpy(float)
            years = 2000 + (jds - 2451545.0) / 365.25
            dt = years - t['COORD_EPOCH']
            pra = t['RA'] + t['PMRA'] * dt / (3.6e6 * np.cos(np.deg2rad(t['DEC'])))
            pdec = t['DEC'] + t['PMDEC'] * dt / 3.6e6
            dcat = sep_arcsec(ra, dec, t['RA'], t['DEC'])
            depoch = sep_arcsec(ra, dec, pra, pdec)
            dist = np.minimum(dcat, depoch)
            candidates = np.where(finite_jd & np.isfinite(dist) & (dist <= threshold))[0]
            if not len(candidates):
                unmatched += 1
                errors.append(f'{filename}: nincs közeli forrás')
                continue
            if not is_multi:
                # Single-exposure files contain many stars: use closest one.
                candidates = [candidates[np.argmin(dist[candidates])]]
            else:
                multi_rows += len(candidates)
            for idx in candidates:
                r = df.iloc[idx]
                def num(col):
                    return pd.to_numeric(r.get(col, np.nan), errors='coerce')
                mag, merr, snr = num('Mag(obs)'), num('Mag(err)'), num('SNR')
                status = 'OK'
                if not np.isfinite(mag):
                    status = 'missing_mag'
                elif not np.isfinite(merr) or merr <= 0:
                    status = 'missing_or_invalid_error'
                elif merr > 0.20:
                    status = 'high_error'
                elif not np.isfinite(snr) or snr < 10:
                    status = 'low_snr'
                rows.append(dict(filename=Path(filename).name, JD=float(jds[idx]),
                    exp_start=str(r.get('EXP start', '')), RA=num('RA'), DEC=num('DEC'),
                    is_variable=str(r.get('Is Variable', '')), mag=mag,
                    mag_err=merr, SNR=snr, reference_stars=num('Valid Reference Stars'),
                    sep_gaia2016_arcsec=float(dcat[idx]), sep_epoch_arcsec=float(depoch[idx]),
                    selected_match='Gaia2016' if dcat[idx] <= depoch[idx] else 'observation_epoch',
                    nearest_sep_arcsec=float(dist[idx]), status=status, valid=status=='OK',
                    source_format='combined' if is_multi else 'individual'))
        except Exception as e:
            errors.append(f'{filename}: {type(e).__name__}: {e}')
    if not rows:
        raise RuntimeError('A célcsillag nem található. Ellenőrizd a változó koordinátáit.')
    data = pd.DataFrame(rows).sort_values('JD').reset_index(drop=True)
    # Remove overlapping exports of the very same exposure. Round only for
    # duplicate grouping; preserve original JDs in all returned measurements.
    # 0.5 second tolerance accommodates rounding, but not 10 s exposures.
    key = np.rint(data.JD.to_numpy() * 86400 * 2).astype('int64')
    data['_jd_key'] = key
    data['_quality_rank'] = (~data.valid).astype(int)
    data['_err_rank'] = data.mag_err.where(data.mag_err > 0, np.inf).fillna(np.inf)
    data = (data.sort_values(['_jd_key', '_quality_rank', '_err_rank'])
                .drop_duplicates('_jd_key', keep='first')
                .drop(columns=['_jd_key', '_quality_rank', '_err_rank'])
                .sort_values('JD').reset_index(drop=True))
    duplicates = len(rows) - len(data)
    print(f'Összefűzött CSV-k: {multi_files}; ezekből koordinátával egyező sorok: {multi_rows}')
    print(f'Ismétlődő időpont miatt kihagyva: {duplicates}; egyedi mérések: {len(data)}')
    return data, seen, unmatched, errors

def add_phase(data,t):
    if t['T0_TYPE']!='HJD':
        raise ValueError(f'Az egyelőre támogatott T0_TYPE: HJD; kapott: {t["T0_TYPE"]}')
    from astropy.time import Time
    from astropy.coordinates import SkyCoord, EarthLocation
    import astropy.units as u
    # Heliocentric correction. Assume both CSV JD and VSX HJD use UTC;
    # the VSX time scale must be verified for rigorous minima timings.
    dt=pd.to_numeric(data.JD,errors='coerce').to_numpy(float)
    if t['LOCATION_MODE'] == 'observatory':
        location=EarthLocation.from_geodetic(lon=t['LONGITUDE']*u.deg,
                                            lat=t['LATITUDE']*u.deg,
                                            height=t['ELEVATION']*u.m)
    else:
        # Fallback: Earth center, with a maximum topocentric difference of ~21 ms.
        location=EarthLocation.from_geocentric(0*u.m,0*u.m,0*u.m)
    time=Time(dt,format='jd',scale='utc',location=location)
    sky=SkyCoord(ra=t['RA']*u.deg,dec=t['DEC']*u.deg,frame='icrs',
                 pm_ra_cosdec=t['PMRA']*u.mas/u.yr,pm_dec=t['PMDEC']*u.mas/u.yr,
                 obstime=Time(t['COORD_EPOCH'],format='jyear',scale='tdb'))
    # Use fixed Gaia coordinates for heliocentric correction; shift from proper
    # motion over 10 years is < 0.002 deg and negligible for light time.
    hjd=time.utc.jd+time.light_travel_time(sky,kind='heliocentric').to_value(u.day)
    # VSX HJD time scale often unspecified, so preserve this assumption explicitly.
    data['HJD_UTC_approx']=hjd
    data['phase']=((hjd-t['T0'])/t['PERIOD'])%1
    data['phase_centered']=((data['phase']+0.5)%1)-0.5
    return data

def observing_sessions(good):
    """Group valid observations into observing nights (UTC noon-to-noon)."""
    utc = pd.to_datetime(good['JD'], unit='D', origin='julian', utc=True)
    # The noon boundary keeps measurements on both sides of midnight
    # in the same observing night, without merging consecutive nights.
    night = (utc - pd.Timedelta(hours=12)).dt.strftime('%Y-%m-%d')
    frames = []
    for label, group in good.assign(_utc=utc, _night=night).groupby('_night', sort=True):
        frames.append((label, group._utc.min(), group._utc.max(), len(group)))
    return frames



def compressed_night_axis(good, separator_hours=0.55):
    """Display observing nights adjacent to each other, without altering JD data."""
    utc = pd.to_datetime(good['JD'], unit='D', origin='julian', utc=True)
    labels = (utc - pd.Timedelta(hours=12)).dt.strftime('%Y-%m-%d')
    jd = good['JD'].to_numpy(dtype=float)
    display = np.zeros(len(jd), dtype=float)
    sessions = []
    cursor = 0.0
    for label in sorted(labels.unique()):
        mask = (labels == label).to_numpy()
        indices = np.where(mask)[0]
        start_jd = jd[indices].min()
        width = (jd[indices].max() - start_jd) * 24
        display[indices] = cursor + (jd[indices] - start_jd) * 24
        sessions.append((label, cursor, cursor + width))
        cursor += width + separator_hours
    return display, sessions


def mark_night_breaks(ax, sessions, offset, jd_values, display_values):
    """Label the compressed observing-night axis with REAL Julian Dates.

    Spacing between nights is deliberately not proportional to elapsed JD.
    """
    tick_positions, tick_labels = [], []
    n = len(sessions)
    per_night = 4 if n <= 3 else (3 if n <= 7 else 2)
    for label, start, end in sessions:
        # The display axis is measured in observing hours, but each tick
        # corresponds to the original JD within its particular night.
        jd_start = float(jd_values[np.argmin(np.abs(display_values - start))])
        for fraction in np.linspace(0, 1, per_night):
            pos = start + fraction * (end - start)
            tick_positions.append(pos)
            tick_labels.append(f'{jd_start + (pos-start)/24 - offset:.3f}')
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, fontsize=8 if n <= 4 else 7,
                       rotation=35 if n > 2 else 0, ha='right' if n > 2 else 'center')
    for (_, _, end), (_, start, _) in zip(sessions, sessions[1:]):
        ax.axvline((end + start) / 2, color='0.70', lw=.8, ls=':', alpha=.6)
    ax.set_xlabel(f'JD − {offset:.1f} UTC')


def flag_spikes(data, window_minutes=20.0, bright_threshold=0.22, sigma=3.0):
    """Conservative candidate filter for brief *brightening* artifacts.

    Estimates a fainter reference track with the 70th percentile over a
    20-minute real-JD window. Unlike v2 it does NOT cancel candidates merely
    because several neighboring points are bad. Strong eclipses (fading)
    are never removed by this asymmetric filter. Keep every original row
    and the local diagnostic values for review.
    """
    data = data.copy()
    data['spike_rejected'] = False
    data['spike_local_median'] = np.nan  # legacy field: reference percentile
    data['spike_deviation_mag'] = np.nan
    data['spike_threshold_mag'] = np.nan
    ix = data.index[data.valid].to_numpy()
    if len(ix) < 15:
        return data
    times = data.loc[ix, 'JD'].to_numpy(float)
    mags = data.loc[ix, 'mag'].to_numpy(float)
    errors = data.loc[ix, 'mag_err'].to_numpy(float)
    radius = window_minutes / (2.0 * 1440.0)
    rejected = np.zeros(len(ix), dtype=bool)
    baseline_values = np.full(len(ix), np.nan)
    deviations = np.full(len(ix), np.nan)
    thresholds = np.full(len(ix), np.nan)
    for k, (jd, mag) in enumerate(zip(times, mags)):
        lo = np.searchsorted(times, jd-radius, side='left')
        hi = np.searchsorted(times, jd+radius, side='right')
        neighbor = np.arange(lo, hi)
        neighbor = neighbor[neighbor != k]
        if len(neighbor) < 10:
            continue
        local = mags[neighbor]
        baseline = float(np.quantile(local, 0.70))
        reference = local[local >= np.quantile(local, 0.40)]
        if len(reference) < 6:
            continue
        scatter = max(0.025, 1.4826*np.median(np.abs(reference - np.median(reference))))
        threshold = max(bright_threshold, sigma*scatter, 3.0*errors[k])
        baseline_values[k] = baseline
        deviations[k] = mag - baseline
        thresholds[k] = threshold
        rejected[k] = baseline - mag > threshold
    data.loc[ix, 'spike_rejected'] = rejected
    data.loc[ix, 'spike_local_median'] = baseline_values
    data.loc[ix, 'spike_deviation_mag'] = deviations
    data.loc[ix, 'spike_threshold_mag'] = thresholds
    return data

def output_plots(data,t,out, time_only=False, show_median=True):
    cleaned = data[data.valid & ~data.spike_rejected].copy()
    # RAW means no photometric quality or spike filtering: every finite
    # magnitude is plotted, regardless of status, SNR or mag_err.
    # Non-numeric/absent magnitudes cannot be displayed on a magnitude axis.
    raw = data[np.isfinite(pd.to_numeric(data['JD'], errors='coerce')) &
               np.isfinite(pd.to_numeric(data['mag'], errors='coerce'))].copy()
    rejected = data[data.valid & data.spike_rejected].copy()
    if cleaned.empty:
        print('Nincs elfogadott magnitúdómérés; csak a táblázat és a napló készül.')
        return

    x = cleaned.JD.to_numpy(float)
    m = cleaned.mag.to_numpy(float)
    err = cleaned.mag_err.to_numpy(float)
    offset = np.floor(x.min()-0.5)+0.5
    site = t.get('SITE','').strip() or 'Not specified'
    nights = observing_sessions(cleaned)
    xplot, night_ranges = compressed_night_axis(cleaned)

    # Rejected overlay for the compressed time axis.
    rej_time = rejected.copy()
    if not rej_time.empty:
        session_map = {label: start for label, start, end in night_ranges}
        clean_labels = (pd.to_datetime(cleaned['JD'], unit='D', origin='julian', utc=True) - pd.Timedelta(hours=12)).dt.strftime('%Y-%m-%d')
        night_min = cleaned.assign(_night=clean_labels).groupby('_night')['JD'].min().to_dict()
        keep = []
        xrej = []
        for idx, row in rej_time.iterrows():
            label = (pd.to_datetime(row['JD'], unit='D', origin='julian', utc=True) - pd.Timedelta(hours=12)).strftime('%Y-%m-%d')
            if label not in session_map or label not in night_min:
                continue
            keep.append(idx)
            xrej.append(session_map[label] + (float(row['JD']) - night_min[label]) * 24.0)
        rej_time = rej_time.loc[keep].copy()
        rej_time['xplot'] = xrej

    # Main time plot: cleaned points + rejected points in light gray.
    height = max(6.2, 3.6 + len(nights)*0.36)
    fig = plt.figure(figsize=(15.2, height))
    ax = fig.add_axes([.065,.12,.62,.78])
    aside = fig.add_axes([.71,.08,.28,.84]); aside.axis('off')
    if not rej_time.empty:
        ax.errorbar(rej_time['xplot'].to_numpy(float), rej_time['mag'].to_numpy(float),
                    yerr=rej_time['mag_err'].to_numpy(float), fmt='none', ecolor='#d4d4d4',
                    alpha=.45, lw=.45, zorder=0)
        ax.scatter(rej_time['xplot'].to_numpy(float), rej_time['mag'].to_numpy(float),
                   s=11, c='#c8c8c8', alpha=.70, zorder=1, label='Rejected spikes')
    ax.errorbar(xplot, m, yerr=err, fmt='none', ecolor='#86b3cf', alpha=.35, lw=.5, zorder=2)
    ax.scatter(xplot, m, s=8, c='#3179a5', alpha=.62, zorder=3, label='Cleaned measurements')
    if show_median and len(x) > 15:
        gaps = np.r_[0, np.cumsum(np.diff(x) > 5/1440)]
        any_smoothing = False
        for gid in np.unique(gaps):
            mask = gaps == gid
            if mask.sum() < 7:
                continue
            yy = pd.Series(m[mask]).rolling(15, center=True, min_periods=5).median().to_numpy()
            ax.plot(xplot[mask], yy, c='#d97915', lw=1.6, alpha=.95,
                    label='Running median' if not any_smoothing else None, zorder=4)
            any_smoothing = True
    ax.invert_yaxis(); ax.set_title(f'{t["NAME"]} — Seestar S50 photometry', fontsize=15)
    mark_night_breaks(ax, night_ranges, offset, x, xplot)
    ax.set_ylabel('Seestar observed magnitude')
    ax.grid(alpha=.16); ax.legend(loc='best')
    object_info = ''
    if time_only and t.get('OBJECT_INFO'):
        obj = t['OBJECT_INFO']
        object_info = ('OBJECT\n' + '\n'.join(f'{key}: {val}' for key,val in obj.items() if val not in ('', None)) + '\n\n')
    info=(object_info + f'OBSERVATION\nSite: {site}\n'
          f'Raw plotted (finite mag): {len(raw)} / {len(data)}\n'
          f'Cleaned measurements: {len(cleaned)}\n'
          f'Rejected spikes shown in gray: {len(rejected)}\n'
          f'Observing nights: {len(nights)}\n'
          f'JD span: {x.max()-x.min():.4f} d\n\n'
          f'PHOTOMETRY\n'
          f'Mag (2–98%): {np.quantile(m,.02):.3f}–{np.quantile(m,.98):.3f}\n'
          f'Median error: {np.median(err):.3f} mag\n'
          f'Median SNR: {np.nanmedian(cleaned.SNR):.1f}\n\n'
          'OBSERVING SESSIONS (UT)\n')
    total_session_hours = sum((end-start).total_seconds()/3600 for _,start,end,_ in nights)
    info += f'Total session span: {total_session_hours:.2f} h\n\n'
    for label,start,end,n in nights:
        start_str=start.strftime('%Y-%m-%d %H:%M')
        end_str=end.strftime('%Y-%m-%d %H:%M')
        elapsed_hours = (end-start).total_seconds()/3600
        info += f'{start_str} –\n  {end_str}  |  {elapsed_hours:.2f} h  |  {n} cleaned points\n'
    font=9.1 if len(nights)<15 else 8.5
    aside.text(0,.99,info,va='top',fontsize=font,linespacing=1.38, transform=aside.transAxes)
    fig.savefig(out/f'{t["NAME"].replace(" ","_")}_time_lightcurve.png', dpi=200, bbox_inches='tight')
    plt.close(fig)
    if time_only:
        return

    # RAW plot: no quality filtering, no spike filtering and no status coloring.
    # Only rows without a numerical magnitude/JD are impossible to display.
    rxplot, raw_nights = compressed_night_axis(raw)
    raw_offset = np.floor(raw['JD'].min()-0.5)+0.5
    fig, ax = plt.subplots(figsize=(11, 5.8))
    # Error bars are shown only where they are finite and positive; the
    # underlying point is ALWAYS drawn even if its error estimate is invalid.
    raw_err = pd.to_numeric(raw['mag_err'], errors='coerce').to_numpy(float)
    has_err = np.isfinite(raw_err) & (raw_err > 0)
    if has_err.any():
        ax.errorbar(rxplot[has_err], raw['mag'].to_numpy(float)[has_err],
                    yerr=raw_err[has_err], fmt='none', ecolor='#a9bdc7',
                    elinewidth=.45, capsize=0, alpha=.60, zorder=1)
    ax.scatter(rxplot, raw['mag'].to_numpy(float), s=9, c='#367aa2',
               alpha=.72, zorder=2)
    ax.invert_yaxis()
    ax.set_title(f'{t["NAME"]} — raw photometry / Seestar S50')
    mark_night_breaks(ax, raw_nights, raw_offset,
                      raw['JD'].to_numpy(float), rxplot)
    ax.set_ylabel('Seestar observed magnitude')
    ax.grid(alpha=.15)
    fig.text(.99,.01,
             f'Observing site: {site}  |  RAW: {len(raw)} finite magnitude measurements (unfiltered)',
             ha='right', va='bottom', fontsize=9)
    fig.tight_layout(rect=[0,.035,1,1])
    fig.savefig(out/f'{t["NAME"].replace(" ","_")}_raw_lightcurve.png', dpi=200)
    plt.close(fig)

    if 'phase_centered' in cleaned.columns:
        fig,ax=plt.subplots(figsize=(9.8,5.8))
        if not rejected.empty and 'phase_centered' in rejected.columns:
            ax.scatter(rejected['phase_centered'], rejected['mag'], s=12, c='#c8c8c8',
                       alpha=.70, zorder=1, label='Rejected spikes')
        pts=ax.scatter(cleaned.phase_centered, cleaned.mag, c=cleaned.mag_err,
                       cmap='viridis_r', s=13, alpha=.75, zorder=2, label='Cleaned measurements')
        ax.invert_yaxis(); ax.set_xlim(-.5,.5)
        ax.set_xlabel(f'Phase  |  T0 = {t["T0"]:.5f} {t["T0_TYPE"]}  |  P = {t["PERIOD"]:.7f} d')
        ax.set_ylabel('Seestar observed magnitude')
        ax.set_title(f'{t["NAME"]} — phase light curve / Seestar S50')
        fig.colorbar(pts,ax=ax,label='Magnitude error')
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend(loc='best')
        ax.grid(alpha=.16)
        fig.text(.99,.01,f'Observing site: {site}',ha='right',va='bottom',fontsize=9)
        fig.tight_layout(rect=[0,.035,1,1])
        fig.savefig(out/f'{t["NAME"].replace(" ","_")}_phase_lightcurve.png', dpi=200)
        plt.close(fig)

def choose_target(valid_configs):
    """Return selected (path, ConfigParser), or None when cancelled."""
    if len(valid_configs) == 1:
        return valid_configs[0]
    try:
        import tkinter as tk
        from tkinter import ttk
        root = tk.Tk()
        root.title('Seestar Photometry – csillagválasztás')
        root.geometry('430x340')
        root.minsize(350, 240)
        root.attributes('-topmost', True)
        frame = ttk.Frame(root, padding=14)
        frame.pack(fill='both', expand=True)
        ttk.Label(frame, text='Válaszd ki a feldolgozandó változót:').pack(anchor='w', pady=(0, 9))
        row = ttk.Frame(frame)
        row.pack(fill='both', expand=True)
        scrollbar = ttk.Scrollbar(row, orient='vertical')
        stars = tk.Listbox(row, height=11, exportselection=False, yscrollcommand=scrollbar.set)
        scrollbar.config(command=stars.yview)
        stars.pack(side='left', fill='both', expand=True)
        scrollbar.pack(side='right', fill='y')
        for f, cp in valid_configs:
            stars.insert(tk.END, cp.sections()[0])
        stars.selection_set(0)
        stars.focus_set()
        selection = [None]
        def confirm(event=None):
            indices = stars.curselection()
            if indices:
                selection[0] = valid_configs[indices[0]]
                root.destroy()
        def cancel(event=None):
            root.destroy()
        buttons = ttk.Frame(frame)
        buttons.pack(fill='x', pady=(12, 0))
        ttk.Button(buttons, text='Feldolgozás', command=confirm).pack(side='right')
        ttk.Button(buttons, text='Mégse', command=cancel).pack(side='right', padx=(0, 8))
        stars.bind('<Double-Button-1>', confirm)
        root.bind('<Return>', confirm)
        root.bind('<Escape>', cancel)
        root.mainloop()
        return selection[0]
    except Exception as exc:
        print(f'Grafikus csillagválasztó nem elérhető: {exc}')
        for i, (_, cp) in enumerate(valid_configs, 1):
            print(f'{i}. {cp.sections()[0]}')
        try:
            answer = input('Add meg a csillag sorszámát (Enter = kilépés): ').strip()
            if not answer:
                return None
            index = int(answer) - 1
            return valid_configs[index] if 0 <= index < len(valid_configs) else None
        except (ValueError, EOFError):
            return None


def main():
    p=argparse.ArgumentParser(description='Seestar CSV light curves')
    p.add_argument('--folder',type=Path)
    p.add_argument('--out',type=Path)
    p.add_argument('--star',type=str,help='Csillagnév kiválasztása ablak nélkül (pl. --star "VW Cep")')
    p.add_argument('--max-sep',type=float,default=8.0,help='Matching radius, arcseconds')
    args=p.parse_args()
    if not args.folder:
        from tkinter import Tk, filedialog
        root=Tk();root.withdraw();root.attributes('-topmost',True)
        chosen=filedialog.askdirectory(title='Válaszd ki a Seestar CSV könyvtárát')
        root.destroy()
        if not chosen:print('Nem választottál könyvtárat.');return
        args.folder=Path(chosen)
    if not args.folder.is_dir():p.error(f'A könyvtár nem létezik: {args.folder}')
    config_files=[f for f in args.folder.glob('*.txt') if f.name.lower() not in ('readme.txt',)]
    valid_configs=[]
    for f in config_files:
        try:
            cp=read_targets(f)
            if len(cp.sections()) == 1 and f.stem.strip().casefold() == cp.sections()[0].strip().casefold():
                valid_configs.append((f,cp))
        except Exception:
            continue
    valid_configs.sort(key=lambda item: item[1].sections()[0].casefold())
    if not valid_configs:
        raise SystemExit('Nem található érvényes csillagnév.txt a kijelölt mappában. A fájl neve és [fejléce] egyezzen.')
    if args.star:
        matching = [item for item in valid_configs if item[1].sections()[0].casefold() == args.star.strip().casefold()]
        if not matching:
            raise SystemExit('Nincs ilyen csillag: '+args.star+'. Választható: '+', '.join(cp.sections()[0] for _, cp in valid_configs))
        selected = matching[0]
    else:
        selected = choose_target(valid_configs)
    if selected is None:
        print('Nem választottál csillagot.')
        return
    cfg,cp=selected
    t=get_target(cp,cp.sections()[0])
    print(f'Célcsillag: {t["NAME"]}; adatlap: {cfg.name}')
    # names=sorted(f for f in args.folder.rglob('*.csv') if 'eredmenyek' not in [x.lower() for x in f.parts])
    names=sorted(args.folder.glob('*.csv'))
    data,n,unmatched,errors=collect(((str(f),f.read_bytes()) for f in names),t,args.max_sep)
    out=args.out or args.folder/'Eredmenyek'
    out.mkdir(parents=True,exist_ok=True)
    phase_message=''
    location_description=(f"{t['SITE'] or 'Unnamed site'}; lat={t['LATITUDE']:.6f} deg, "
                          f"lon={t['LONGITUDE']:.6f} deg, elevation={t['ELEVATION']:.1f} m"
                          if t['LOCATION_MODE']=='observatory' else 'Geocenter (no site coordinates in TXT)')
    print(f'Heliocentric correction location: {location_description}')
    try:data=add_phase(data,t)
    except ImportError:
        phase_message='Fázisgörbe nem készül: telepítsd az astropy csomagot (python -m pip install astropy).'
    except Exception as e:phase_message=f'Fázisgörbe nem készült: {e}'
    base=t['NAME'].replace(' ','_')
    data=flag_spikes(data)
    data.to_csv(out/f'{base}_all_measurements.csv',index=False,float_format='%.9f')
    data[data.valid & ~data.spike_rejected].to_csv(out/f'{base}_cleaned.csv',index=False,float_format='%.9f')
    data[data.spike_rejected].to_csv(out/f'{base}_rejected_spikes.csv',index=False,float_format='%.9f')
    print(f'Automatikusan kiszűrt kiugró pontok: {int(data.spike_rejected.sum())}; tisztított mérések: {int((data.valid & ~data.spike_rejected).sum())}')
    output_plots(data,t,out)
    counts=data.status.value_counts().to_string()
    report=(f'Object: {t["NAME"]}\nCSV files: {n}\nMatched files: {len(data)}\nUnmatched files: {unmatched}\n'
            f'Other processing errors: {len(errors)}\nValid magnitudes (raw): {int(data.valid.sum())}\n'
            f'Rejected brightening candidates: {int(data.spike_rejected.sum())}\nCleaned magnitudes: {int((data.valid & ~data.spike_rejected).sum())}\n'
            f'Identification: nearest of Gaia DR3 J2016 position and proper-motion-corrected position; maximum {args.max_sep}".\n'
            f'Gaia2016 matches: {sum(data.selected_match=="Gaia2016")}\n'
            f'Observation-epoch matches: {sum(data.selected_match=="observation_epoch")}\n'
            f'\nStatus counts:\n{counts}\n\n{phase_message}\n'
            'CAVEAT: Gaia/CSV coordinate epoch may differ. HJD correction assumes CSV JD=UTC; VSX HJD time scale should be verified.\n'
            '\nFirst processing messages:\n'+'\n'.join(errors[:100]))
    (out/f'{base}_quality_report.txt').write_text(report,encoding='utf-8')
    print(report.split('First processing messages:')[0]);print(f'Eredmények: {out}')

if __name__=='__main__':main()

# ---------------------------------------------------------------------------
# All Variables mode: time plots only, no phase/RAW plots or output tables.
# ---------------------------------------------------------------------------
def all_variables(folder, match_arcsec=2.0, vsx_radius_arcsec=5.0,
                  min_points=100, show_median=True):
    """Generate one time PNG for each coordinate-matched Seestar variable.

    The input directory is intentionally non-recursive. VSX names are optional:
    when offline, not found or ambiguous, a coordinate-based identifier is used.
    """
    import re
    from scipy.spatial import cKDTree

    if int(min_points) < 2:
        raise ValueError('A minimális pontszám legalább 2 legyen.')
    min_points = int(min_points)
    folder = Path(folder)
    files = sorted(folder.glob('*.csv'))
    if not files:
        raise ValueError('Nincs CSV-fájl a kiválasztott könyvtárban.')
    output = folder / 'Eredmenyek' / 'All_Variables'
    output.mkdir(parents=True, exist_ok=True)
    observations, file_errors = [], []
    required = {'RA', 'DEC', 'Is Variable', 'Mag(obs)', 'Mag(err)', 'EXP start(JD)'}
    print(f'All Variables: {len(files)} CSV fájl, almappák nélkül.')
    for i, file in enumerate(files, 1):
        try:
            table = pd.read_csv(file, usecols=lambda c: c in required or c == 'SNR', low_memory=False)
            if not required.issubset(table.columns):
                file_errors.append(f'{file.name}: hiányzó oszlopok')
                continue
            sel = table['Is Variable'].astype(str).str.strip().str.lower().isin(['true', '1', 'yes'])
            part = table.loc[sel].copy()
            for key in ['RA','DEC','Mag(obs)','Mag(err)','EXP start(JD)','SNR']:
                if key in part:
                    part[key] = pd.to_numeric(part[key], errors='coerce')
            part = part.rename(columns={'Mag(obs)':'mag', 'Mag(err)':'mag_err', 'EXP start(JD)':'JD'})
            part = part.dropna(subset=['RA','DEC','JD'])
            part = part[(part.DEC.between(-90, 90)) & (part.RA.between(0,360, inclusive='left'))]
            observations.append(part[['RA','DEC','JD','mag','mag_err'] + (['SNR'] if 'SNR' in part.columns else [])])
        except Exception as exc:
            file_errors.append(f'{file.name}: {exc}')
        if i % 250 == 0:
            print(f'Beolvasva: {i}/{len(files)} fájl')
    if not observations or not any(len(d) for d in observations):
        raise RuntimeError('Az Is Variable=True jelölésű csillagokhoz nem találtam méréseket.')
    df = pd.concat(observations, ignore_index=True)
    print(f'Változójelölésű sorok: {len(df)}, hibás fájlok: {len(file_errors)}')

    ra = np.deg2rad(df.RA.to_numpy(float))
    dec = np.deg2rad(df.DEC.to_numpy(float))
    xyz = np.stack([np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)],axis=1)
    # The catalog has repeated measurements, so group connected points using
    # an angular nearest-neighbour graph. Record ambiguous matches separately.
    tree = cKDTree(xyz)
    close = tree.query_pairs(2*np.sin(np.deg2rad(match_arcsec/3600)/2), output_type='ndarray')
    parent = np.arange(len(df))
    rank = np.zeros(len(df), dtype=np.uint8)
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(a,b):
        a,b=find(a),find(b)
        if a==b:return
        if rank[a]<rank[b]:a,b=b,a
        parent[b]=a
        if rank[a]==rank[b]:rank[a]+=1
    for a,b in close:union(int(a),int(b))
    groups = {}
    for i in range(len(df)):
        root=find(i)
        groups.setdefault(root, []).append(i)
    groups = sorted(groups.values(), key=len, reverse=True)
    print(f'Koordináta szerint különböző jelöltek: {len(groups)}')

    # The VSX api.list radius parameter is in DEGREES, not arcseconds.
    # Grouping (match_arcsec) and catalog identification (vsx_radius_arcsec)
    # intentionally remain independent.
    from urllib.parse import urlencode
    from urllib.request import Request, urlopen
    import json
    from urllib.error import URLError

    def parse_vsx(body, ra0, dec0):
        payload = json.loads(body)
        obj = payload.get('VSXObjects', payload)
        if isinstance(obj, dict):
            objects = obj.get('VSXObject', [])
        elif isinstance(obj, list):
            objects = obj
        else:
            objects = []
        if isinstance(objects, dict):
            objects = [objects]
        matches = []
        for item in objects:
            if not isinstance(item, dict):
                continue
            try:
                vra = float(item.get('RA2000'))
                vdec = float(item.get('Declination2000'))
                separation = float(sep_arcsec(vra, vdec, ra0, dec0))
                if separation <= vsx_radius_arcsec + 0.00001:
                    name = str(item.get('Name', '')).strip()
                    if name:
                        period = item.get('Period', '')
                        vartype = item.get('Type', '')
                        matches.append(dict(sep=separation, name=name, ra=vra,
                                            dec=vdec, period=period, vartype=vartype))
            except (ValueError, TypeError):
                continue
        return sorted(matches, key=lambda x:x['sep'])

    def format_coords(ra, dec):
        rah = (ra % 360) / 15
        hh = int(rah); mm = int((rah-hh)*60)
        ss = ((rah-hh)*60-mm)*60
        d = abs(dec); dd = int(d); dm = int((d-dd)*60)
        ds = ((d-dd)*60-dm)*60
        return f'{hh:02d}:{mm:02d}:{ss:05.2f}', f'{"+" if dec >= 0 else "-"}{dd:02d}:{dm:02d}:{ds:04.1f}'

    def vsx_lookup(ra0, dec0):
        params = urlencode({'view': 'api.list', 'ra': f'{ra0:.8f}',
                            'dec': f'{dec0:.8f}',
                            'radius': f'{vsx_radius_arcsec / 3600:.9f}',
                            'format': 'json'})
        url = 'https://www.aavso.org/vsx/index.php?' + params
        req = Request(url, headers={'User-Agent': 'SeestarPhotometry/1.1 (VSX coordinate query)',
                                    'Accept': 'application/json'})
        with urlopen(req, timeout=15) as response:
            body = response.read().decode('utf-8-sig')
        return parse_vsx(body, ra0, dec0)

    print(f'VSX koordináta-keresés: {vsx_radius_arcsec:g}″; csillagcsoportosítás: {match_arcsec:g}″')
    vsx_available = True
    vsx_matched = vsx_missing = vsx_ambiguous = 0

    import matplotlib.pyplot as plt
    made = 0
    used_names = set()
    summary = []
    skipped = 0
    for n, ix in enumerate(groups, 1):
        d = df.iloc[ix].copy().sort_values('JD')
        rvec = np.mean(xyz[ix], axis=0)
        rvec /= np.linalg.norm(rvec)
        cra = np.rad2deg(np.arctan2(rvec[1],rvec[0])) % 360
        cdec = np.rad2deg(np.arcsin(rvec[2]))
        fallback = f'VAR_RA{cra:09.5f}_DEC{cdec:+09.5f}'
        # Filter on unique, valid photometric samples BEFORE VSX network calls.
        finite = np.isfinite(d.mag.to_numpy(float)) & np.isfinite(d.JD.to_numpy(float))
        d = d.loc[finite].copy()
        if d.empty:
            skipped += 1
            continue
        d['_time_key'] = np.rint(d.JD.to_numpy(float)*86400*2).astype('int64')
        d['_err_sort'] = d.mag_err.where(d.mag_err.gt(0), np.inf).fillna(np.inf)
        d = (d.sort_values(['_time_key','_err_sort']).drop_duplicates('_time_key')
               .drop(columns=['_time_key','_err_sort']).sort_values('JD'))
        err_pre = pd.to_numeric(d.mag_err, errors='coerce').to_numpy(float)
        snr_pre = pd.to_numeric(d.SNR, errors='coerce').to_numpy(float) if 'SNR' in d else np.full(len(d), np.nan)
        pre_valid = np.isfinite(err_pre) & (err_pre>0) & (err_pre<=0.2) & np.isfinite(snr_pre) & (snr_pre>=10)
        if int(pre_valid.sum()) < min_points:
            skipped += 1
            continue
        name = fallback
        vsx_obj = None
        vsx_status = 'not_found'
        if vsx_available:
            try:
                matches = vsx_lookup(cra, cdec)
                if len(matches) == 1:
                    vsx_obj = matches[0]
                    name = vsx_obj['name']
                    vsx_status = 'matched'
                    vsx_matched += 1
                    print(f"VSX: {fallback} -> {name} ({vsx_obj['sep']:.2f}″)")
                elif len(matches) > 1:
                    vsx_ambiguous += 1
                    vsx_status = 'ambiguous'
                    print(f'VSX több találat ({fallback}): ' + ', '.join(x['name'] for x in matches[:6]) + '; koordinátanév marad.')
                else:
                    vsx_missing += 1
                    print(f'VSX nincs találat: {fallback}')
            except Exception as exc:
                print(f'VSX lekérdezés sikertelen ({fallback}): {type(exc).__name__}: {str(exc)[:160]}')
                print('A következő csillagok koordinátanévvel készülnek. Ellenőrizd az internetkapcsolatot és a VSX API-t.')
                vsx_available = False
                vsx_status = 'query_failed'
        err = pd.to_numeric(d.mag_err, errors='coerce').to_numpy(float)
        snr = pd.to_numeric(d.SNR, errors='coerce').to_numpy(float) if 'SNR' in d else np.full(len(d), np.nan)
        valid = np.isfinite(err) & (err>0) & (err<=0.2) & np.isfinite(snr) & (snr>=10)
        # Reuse exactly the target star's time plotting engine, without RAW/phase plots.
        clean_name = re.sub(r'[^a-zA-Z0-9_+.-]+', '_', name).strip('._')[:100] or fallback
        if clean_name in used_names:
            clean_name += f'_{n}'
        used_names.add(clean_name)
        d = d.reset_index(drop=True)
        d['SNR'] = snr
        d['status'] = np.select(
            [~np.isfinite(d.mag.to_numpy(float)),
             ~np.isfinite(err) | (err <= 0),
             err > 0.2,
             ~np.isfinite(snr) | (snr < 10)],
            ['missing_mag', 'missing_or_invalid_error', 'high_error', 'low_snr'],
            default='OK')
        d['valid'] = d['status'].eq('OK')
        d = flag_spikes(d)
        good = d[d['valid'] & ~d['spike_rejected']]
        if len(good) < min_points:
            skipped += 1
            print(f'{name}: {len(good)} tisztított pont < {min_points}, kihagyva.')
            continue
        coord_ra, coord_dec = (vsx_obj['ra'], vsx_obj['dec']) if vsx_obj else (cra, cdec)
        ra_text, dec_text = format_coords(coord_ra, coord_dec)
        period = vsx_obj.get('period') if vsx_obj else None
        try:
            period_num = float(str(period).replace(',', '.'))
            period_text = f'{period_num:.9g} d ({period_num * 24:.4f} h)' if np.isfinite(period_num) and period_num > 0 else '--'
        except (TypeError, ValueError):
            period_text = '--'
        vartype = vsx_obj.get('vartype') if vsx_obj else None
        object_info = {'Name': name, 'RA (J2000)': ra_text, 'Dec (J2000)': dec_text,
                       'Period': period_text, 'Type': vartype or '--'}
        output_plots(d, {'NAME': clean_name, 'SITE': 'Esztergom, HU',
                         'OBJECT_INFO': object_info}, output, time_only=True,
                     show_median=show_median)
        made += 1
        first = pd.to_datetime(good.JD.min(), unit='D', origin='julian', utc=True)
        last = pd.to_datetime(good.JD.max(), unit='D', origin='julian', utc=True)
        summary.append({'ID': name, 'VSX_status': vsx_status,
                        'RA_deg': round(coord_ra, 7), 'DEC_deg': round(coord_dec, 7),
                        'RA_J2000': ra_text, 'DEC_J2000': dec_text,
                        'Period_d': period if period_text != '--' else '',
                        'Type': vartype or '', 'Points_raw': len(d),
                        'Points_cleaned': len(good), 'Rejected_spikes': int(d.spike_rejected.sum()),
                        'JD_start': float(good.JD.min()), 'JD_end': float(good.JD.max()),
                        'UTC_start': first.strftime('%Y-%m-%d %H:%M:%S'),
                        'UTC_end': last.strftime('%Y-%m-%d %H:%M:%S'),
                        'Span_hours': round((good.JD.max()-good.JD.min())*24, 4),
                        'Median_mag': round(float(good.mag.median()), 4),
                        'Median_error_mag': round(float(good.mag_err.median()), 4),
                        'Median_SNR': round(float(good.SNR.median()), 2),
                        'Plot_file': f'{clean_name}_time_lightcurve.png'})
        print(f'{n}/{len(groups)}  {name}: {len(d)} pont')
    columns = ['ID','VSX_status','RA_deg','DEC_deg','RA_J2000','DEC_J2000',
               'Period_d','Type','Points_raw','Points_cleaned','Rejected_spikes',
               'JD_start','JD_end','UTC_start','UTC_end','Span_hours',
               'Median_mag','Median_error_mag','Median_SNR','Plot_file']
    summary_path = output / 'All_Variables_summary.csv'
    pd.DataFrame(summary, columns=columns).to_csv(summary_path, index=False, encoding='utf-8-sig')
    print(f'Összesítés: {summary_path}; minimum: {min_points} pont; kihagyva: {skipped}')
    print(f'VSX: azonosított={vsx_matched}, nincs találat={vsx_missing}, több találat={vsx_ambiguous}.')
    print(f'Elkészült time-görbék: {made}. Mappa: {output}')
    if file_errors:
        print('Első beolvasási hibák:', *file_errors[:10], sep='\n  ')
    return made
