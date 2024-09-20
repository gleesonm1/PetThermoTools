import urllib.request
import zipfile
import os
import sys
import site
import sysconfig
import importlib


def install_alphaMELTS(chip="Linux", file_location = None):
    '''
    Download, extract, and add the alphaMELTS for Python files to the Python path.
    Either store the files in the current working directory or (using the file_location kwarg) add them to a location of your choice.
    chip - string variable with 4 options based on OS and chip type.
        'Apple' - M1, M2, M3 etc. chips for MacOS
        'Intel4Mac' - older MacBooks etc. that use an Intel chip.
        'Windows' - windows operating systems
        'Linux' - Default.
    '''
    try:
        from meltsdynamic import meltsdynamic
        print('alphaMELTS already installed and added to Python path')
        return
    except:
        # URL of the file
        if chip == "Apple":
            url = "https://github.com/magmasource/alphaMELTS/releases/download/v2.3.1/alphamelts-py-2.3.1-macos-arm64.zip"
        elif chip == "Intel4Mac":
            url = "https://github.com/magmasource/alphaMELTS/releases/download/v2.3.1/alphamelts-py-2.3.1-macos-x86_64.zip"
        elif chip == "Windows":
            url = "https://github.com/magmasource/alphaMELTS/releases/download/v2.3.1/alphamelts-py-2.3.1-win64.zip"
        else:
            url = "https://github.com/magmasource/alphaMELTS/releases/download/v2.3.1/alphamelts-py-2.3.1-linux.zip"

        # Path to save the file
        if file_location is None:
            if chip == "Apple":
                zip_path = "alphamelts-py-2.3.1-macos-arm64.zip"
            elif chip == "Intel4Mac":
                zip_path = "alphamelts-py-2.3.1-macos-x86_64.zip"
            elif chip == "Windows":
                zip_path = "alphamelts-py-2.3.1-win64.zip"
            else:
                zip_path = "alphamelts-py-2.3.1-linux.zip"
        else:
            if chip == "Apple":
                zip_path = file_location + "alphamelts-py-2.3.1-macos-arm64.zip"
            elif chip == "Intel4Mac":
                zip_path = file_location + "alphamelts-py-2.3.1-macos-x86_64.zip"
            elif chip == "Windows":
                zip_path = file_location + "alphamelts-py-2.3.1-win64.zip"
            else:
                zip_path = file_location + "alphamelts-py-2.3.1-linux.zip"


        # Download the file with error handling
        try:
            urllib.request.urlretrieve(url, zip_path)
        except Exception as e:
            print(f"Error downloading file: {e}")
            sys.exit(1)
            return

        # Path to extract the contents
        if file_location is None:
            extract_path = "alphamelts_py"
        else:
            extract_path = file_location + "alphamelts_py"

        # Extract the zip file with error handling
        try:
            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                zip_ref.extractall(extract_path)
        except zipfile.BadZipFile as e:
            print(f"Error extracting zip file: {e}")
            sys.exit(1)
            return
        
        # Add the extracted directory to the Python path
        if file_location is None:
            sys.path.append(os.path.join(extract_path, zip_path[:-4]))
        else:
            if chip == "Apple":
                sys.path.append(os.path.join(extract_path,"alphamelts-py-2.3.1-macos-arm64"))
            elif chip == "Intel4Mac":
                sys.path.append(os.path.join(extract_path,"alphamelts-py-2.3.1-macos-x86_64"))
            elif chip == "Windows":
                sys.path.append(os.path.join(extract_path,"alphamelts-py-2.3.1-win64"))
            else:
                sys.path.append(os.path.join(extract_path,"alphamelts-py-2.3.1-linux"))

        try:
            import meltsdynamic
            importlib.reload(meltsdynamic)
            from meltsdynamic import MELTSdynamic
            print('Download and Extraction of alphaMELTS for Python files is successful.')
        except:
            print('Error: alphaMELTS for Python not installed correctly.')
            return
        
        # get current working directory
        cwd = os.getcwd()

        # Get site-packages path
        site_packages_path = sysconfig.get_paths()["purelib"]

        # Define pth file path and custom path
        pth_file_path = os.path.join(site_packages_path, "my_MELTS_path.pth")
        if file_location is None:
            custom_path = os.path.join(cwd, extract_path, zip_path[:-4])
        else:
            custom_path = os.path.join(zip_path[:-4])

        # Normalize paths before writing
        custom_path = os.path.normpath(custom_path)

        # Write to the .pth file with error handling
        try:
            with open(pth_file_path, "w") as f:
                f.write(custom_path)
        except Exception as e:
            print(f"Error writing to .pth file: {e}")
            sys.exit(1)

        return
    
def remove_alphaMELTS_path():
    print('Please note this does not remove the previously downloaded alphaMELTS files. That has to be done normally. \n This function simply removes the Python path to those files so that you can install/update to a new version of alphaMELTS for Python.')
    site_packages_dirs = site.getsitepackages()

    for dir in site_packages_dirs:
        print(f"Checking {dir} for .pth files")
        for file in os.listdir(dir):
            if file.endswith(".pth"):
                print(f"Found .pth file: {file}")

    # Path to your site-packages directory (adjust as needed)
    site_packages_path = site.getsitepackages()[0]

    # Path to the .pth file
    pth_file_path = os.path.join(site_packages_path, "my_MELTS_path.pth")

    # Remove the .pth file
    if os.path.exists(pth_file_path):
        os.remove(pth_file_path)
        print(f"Removed {pth_file_path}")
    else:
        print(f".pth file not found: {pth_file_path}")

    return