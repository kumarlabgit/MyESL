rm -rf _internal
rm MyESL_model_apply.exe
rm MyESL.exe
pyinstaller MyESL.py
pyinstaller MyESL_model_apply.py
mv dist\MyESL\_internal .
mv dist\MyESL\MyESL.exe .
mv dist\MyESL_model_apply\MyESL_model_apply.exe .
rm -rf build
rm -rf dist
rm MyESL_model_apply.spec
rm MyESL.spec