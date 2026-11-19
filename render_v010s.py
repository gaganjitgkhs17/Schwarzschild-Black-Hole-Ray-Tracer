"""Version 0.10S — independent NumPy engine, Matplotlib only for plotting.

Run: python render_v010s.py --size 48 --samples 2
Run local compatibility gate FIRST: python physics_audit_v010s.py
"""
import argparse
import json
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import bh_engine_v010s as engine


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--size',type=int,default=48)
    p.add_argument('--samples',type=int,choices=(1,2,4,8),default=2)
    p.add_argument('--tile-rows',type=int,default=8)
    p.add_argument('--step',type=float,default=.008)
    p.add_argument('--output',type=Path,default=Path('figures'))
    args=p.parse_args()
    rgb,hits=engine.render_linear(args.size,args.samples,args.tile_rows,args.step)
    display=engine.tone_map(rgb)
    args.output.mkdir(parents=True,exist_ok=True)
    stem=args.output/f'disk_standalone_{args.size}px_{args.samples}x{args.samples}_h{args.step:g}'
    np.savez_compressed(str(stem)+'.npz',linear_rgb=rgb,display_rgb=display,hit_fraction=hits)
    fig,ax=plt.subplots(figsize=(8,8))
    ax.imshow(display,origin='lower',extent=[-engine.HALF_WIDTH,engine.HALF_WIDTH]*2,interpolation='nearest')
    ax.set(xlabel='Camera x (M)',ylabel='Camera y (M)',title='Schwarzschild disk — standalone snapshot engine')
    fig.savefig(str(stem)+'.png',dpi=200,bbox_inches='tight')
    plt.close(fig)
    integral=engine.camera_brightness_integral(rgb)
    metadata=dict(version='0.10S',engine_source='archived 0.9A/0.9B snapshot',
                  validated_against_local_legacy=False,
                  brightness_integral=integral,mean_hit_fraction=float(hits.mean()),
                  size=args.size,samples=args.samples,step=args.step,
                  timestamp_utc=datetime.now(timezone.utc).isoformat())
    Path(str(stem)+'.json').write_text(json.dumps(metadata,indent=2),encoding='utf-8')
    print(f'Brightness integral: {integral:.10f}')
    print(f'Subray disk-hit fraction: {hits.mean():.8f}')
    print(f'Saved {stem}.png/.npz/.json')
    print('IMPORTANT: Run physics_audit_v010s.py against YOUR corrected local engine before adoption.')

if __name__=='__main__':main()
