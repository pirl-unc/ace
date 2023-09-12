# -*- mode: python ; coding: utf-8 -*-
from PyInstaller.utils.hooks import copy_metadata

datas = [('/Users/leework/miniconda3/envs/ace/lib/python3.10/site-packages/eel/eel.js', 'eel'), ('views', 'views'), ('trained_model_w_data_augmentation_b3000.pt', '.')]
datas += copy_metadata('tqdm')
datas += copy_metadata('regex')
datas += copy_metadata('filelock')
datas += copy_metadata('requests')
datas += copy_metadata('packaging')
datas += copy_metadata('numpy')
datas += copy_metadata('torch')


block_cipher = None


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
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)
pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='ACE',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=['ace_gui.ico'],
)
coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='ACE',
)
app = BUNDLE(
    coll,
    name='ACE.app',
    icon='ace_gui.ico',
    bundle_identifier=None,
)
