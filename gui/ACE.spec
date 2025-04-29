# -*- mode: python ; coding: utf-8 -*-
from PyInstaller.utils.hooks import copy_metadata

datas = [('/Users/leework/miniconda3/anaconda3/envs/ace/lib/python3.10/site-packages/eel/eel.js', 'eel'), ('views', 'views'), ('trained_model_w_data_augmentation_b3000.pt', '.')]
datas += copy_metadata('tqdm')
datas += copy_metadata('regex')
datas += copy_metadata('filelock')
datas += copy_metadata('requests')
datas += copy_metadata('packaging')
datas += copy_metadata('numpy')
datas += copy_metadata('torch')


a = Analysis(
    ['ace_gui.py'],
    pathex=[],
    binaries=[],
    datas=datas,
    hiddenimports=['bottle_websocket'],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='ACE',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=['ace_gui.ico'],
)
app = BUNDLE(
    exe,
    name='ACE.app',
    icon='ace_gui.ico',
    bundle_identifier=None,
)
