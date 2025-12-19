import urllib.request
import zipfile
import os
import sys
import site
import sysconfig
import importlib
from pathlib import Path
import platform
from petthermotools.Path_wrappers import *
from petthermotools.Compositions import *

def install_MAGEMinCalc():
    '''
    Establish a new julia environment that will be used for any MAGEMin calculations performed through PetThermoTools.
    MAGEMinCalc and MAGEMin_C are added to this new environment. No other installation steps are required.
    '''

    print("If you encounter issues with the installation process please contact me at gleesonm@berkeley.edu")

    from juliacall import Main as jl
    env_dir = Path.home() / ".petthermotools_julia_env"
    env_dir.mkdir(exist_ok=True)
    jl_env_path = env_dir.as_posix()

    jl.seval(f"""
        import Pkg
        Pkg.activate("{jl_env_path}")
        Pkg.instantiate()
        """)
    jl.seval(f"""
        for pkg in values(Pkg.dependencies())
            try
                if pkg.is_direct_dep
                    Pkg.rm(pkg.name)
                end
            catch e
                @warn "Failed to remove package" pkg=pkg.name exception=e
            end
        end
        Pkg.gc()  # Clean up unused dependencies

        """)
    jl.seval(f""" 
        Pkg.Registry.update()
        Pkg.add("PythonCall")
             
        using .Sys
        if Sys.iswindows()
            ENV["JULIA_SSL_CA_ROOTS_PATH"] = ""
        end
             
        # Install MAGEMinCalc and dependencies
        try
            Pkg.add(url="https://github.com/gleesonm1/MAGEMinCalc.git", rev="v0.5.2")
        catch e
            @warn "Failed to install MAGEMinCalc via HTTPS, retrying..." exception=e
            ENV["JULIA_SSL_CA_ROOTS_PATH"] = ""
            Pkg.add(url="https://github.com/gleesonm1/MAGEMinCalc.git", rev="v0.5.2")
        end
             
        Pkg.add(name = "MAGEMin_C", version="2.0.6")
             
        Pkg.resolve()   
        Pkg.precompile()
             
        println("MAGEMin environment ready at {jl_env_path}")
        """)

def update_MAGEMinCalc():
    '''
    Update MAGEMinCalc and MAGEMin_C to the latest versions that are compatible with PetThermoTools. It is useful to run this function once following each upgrade of PetThermoTools to ensure compatibility.
    '''
    from juliacall import Main as jl
    env_dir = Path.home() / ".petthermotools_julia_env"
    jl_env_path = env_dir.as_posix()

    jl.seval(f"""
        import Pkg
        Pkg.activate("{jl_env_path}")

        # Install MAGEMinCalc and dependencies
        try
            Pkg.add(url="https://github.com/gleesonm1/MAGEMinCalc.git", rev="v0.5.2")
        catch e
            @warn "Failed to install MAGEMinCalc via HTTPS, retrying..." exception=e
            ENV["JULIA_SSL_CA_ROOTS_PATH"] = ""
            Pkg.add(url="https://github.com/gleesonm1/MAGEMinCalc.git", rev="v0.5.2")
        end
             
        Pkg.add(name = "MAGEMin_C", version="2.0.6")
             
        Pkg.resolve()   
        Pkg.precompile()
        println("MAGEMin environment updated at {jl_env_path}")
    """)

def test_MAGEMinCalc():
    '''
    Test the MAGEMin installation worked! This function should perform a simple fractional crystallization calculation and print the results.
    '''
    from juliacall import Main as jl
    env_dir = Path.home() / ".petthermotools_julia_env"
    jl_env_path = env_dir.as_posix()

    jl.seval(f"""
        import Pkg
        Pkg.activate("{jl_env_path}")
             
        using MAGEMinCalc

        comp = Dict("SiO2_Liq" => 47.5, "Al2O3_Liq" => 16.4, "CaO_Liq" => 11.6, "MgO_Liq" => 9.38,
                    "FeOt_Liq" => 9.16, "K2O_Liq" => 0.329, "Na2O_Liq" => 2.25, "TiO2_Liq" => 2.29, 
                    "Fe3Fet_Liq" => 0.15, "Cr2O3_Liq" => 0.0, "H2O_Liq" => 0.68)

        Results = MAGEMinCalc.path(comp = comp, T_end_C = 1100.0, dt_C = 2.0, 
                    P_bar = 1000.0, frac_xtal = true, 
                    Model = "ig",
                    find_liquidus = true)

        println(Results["liq1"])

    """)

def install_alphaMELTS(file_location=None, admin=True, version="2.3.1"):
    try:
        from meltsdynamic import meltsdynamic
        print(f'alphaMELTS already installed and available.')
        return
    except ImportError:
        pass

    system = platform.system()
    machine = platform.machine().lower()
    is_arm = any(arch in machine for arch in ["arm", "aarch"])
    
    # 1. Construct the download URL
    if version == "2.3.2":
        # Note: Using the specific naming convention for 2.3.2 assets
        base_url = f"https://github.com/magmasource/alphaMELTS/releases/download/v{version}-beta.0"
        if system == "Darwin":
            tag = "macosx_14_0-arm64" if is_arm else "macosx_14_0-x86_64"
        elif system == "Windows":
            tag = "win64-aarch64" if is_arm else "win64-x86_64"
        elif system == "Linux":
            tag = "ubuntu_22_04-aarch64" if is_arm else "ubuntu_22_04-x86_64"
        else:
            raise OSError(f"Unsupported system: {system}")
        url = f"{base_url}/alphamelts-py-{version}-{tag}.zip"
    elif version == "2.3.1":
        # Fallback for 2.3.1 or older naming conventions
        base_url = f"https://github.com/magmasource/alphaMELTS/releases/download/v{version}"
        if system == "Darwin":
            tag = "macos-arm64" if is_arm else "macos-x86_64"
        elif system == "Windows":
            tag = "win64"
        else:
            tag = "linux"
        url = f"{base_url}/alphamelts-py-{version}-{tag}.zip"
    else:
        raise Warning(f"Version provided not currently supported")

    # 2. Derive paths dynamically from the URL
    zip_name = url.split('/')[-1]
    folder_name = zip_name.replace('.zip', '')
    
    # Ensure file_location is handled safely
    base_dir = file_location if file_location else os.getcwd()
    zip_path = os.path.join(base_dir, zip_name)
    extract_path = os.path.join(base_dir, f"alphamelts_py-{version}")
    full_module_path = os.path.join(extract_path, folder_name)

    # 3. Download
    print(f"Downloading alphaMELTS from: {url}...")
    try:
        urllib.request.urlretrieve(url, zip_path)
    except Exception as e:
        print(f"Error downloading: {e}. Contact Matt Gleeson (gleesonm@berkeley.edu)")
        return

    # 4. Extract
    print(f"Extracting to: {extract_path}...")
    try:
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(extract_path)
    except Exception as e:
        print(f"Extraction failed: {e}")
        return

    # 5. Path Configuration
    if admin:
        # Attempt to add to site-packages via .pth file
        try:
            site_packages = sysconfig.get_paths()["purelib"]
            pth_file = os.path.join(site_packages, "my_MELTS_path.pth")
            with open(pth_file, "w") as f:
                f.write(full_module_path)
            print(f"Path permanently added via: {pth_file}")
            sys.path.append(full_module_path)
        except Exception as e:
            print(f"Admin pathing failed (likely a permissions issue): {e}")
            admin = False # Fall back to printing manual instructions

    if not admin:
        print(f"\nManual Path Configuration Required:")
        print(f"Add these lines to your notebook:")
        print(f"import sys")
        print(f"sys.path.append(r'{full_module_path}')")

    return

# def install_alphaMELTS(file_location = None, admin = True, version = "2.3.1"):
#     '''
#     Downloads, extracts, and optionally adds the platform-specific alphaMELTS for Python files
#     (meltsdynamic.py) to the Python path.

#     Args:
#         file_location (str, optional): Directory to download and extract the files.
#             If None (default), uses the current working directory.
#         admin (bool, optional): If True (default), attempts to permanently add the files to the 
#             Python path using a .pth file (requires admin/root privileges).
#             If False, you must manually append the path in your code (at the start of each notebook using petthermotools).
#         version (str, optional): The alphaMELTS version to download (e.g., "2.3.1").
#             Default is "2.3.1".
    
#     Returns:
#         None. Prints status messages regarding download, extraction, and path configuration.
#     '''
#     try:
#         from meltsdynamic import meltsdynamic
#         print('alphaMELTS already installed and added to Python path')
#         return
#     except:
#         system = platform.system()
#         machine = platform.machine().lower()
#         # version = "2.3.1"

#         if version == "2.3.2":
#             base_url = f"https://github.com/magmasource/alphaMELTS/releases/download/v{version}-beta.0"
#             # macOS
#             if system == "Darwin":
#                 if "arm" in machine:  # Apple Silicon
#                     url = f"{base_url}/alphamelts-py-{version}-macosx_14_0-arm64.zip"
#                 else:  # Intel Mac
#                     url = f"{base_url}/alphamelts-py-{version}-macosx_14_0-x86_64.zip"
#             elif system == "Windows":
#                 if "aarch" in machine:
#                     url = f"{base_url}/alphamelts-py-{version}-win64-aarch64.zip"
#                 else:
#                     url = f"{base_url}/alphamelts-py-{version}-win64-x86_64.zip"
#             elif system == "Linux":
#                 if "aarch" in machine:
#                     url = f"{base_url}/alphamelts-py-{version}-ubuntu_22_04-aarch64.zip"
#                 else:
#                     url = f"{base_url}/alphamelts-py-{version}-ubuntu_22_04-x86_64.zip"
#             else:
#                 raise OSError(f"Unsupported system: {system}")
#         else:
#             base_url = f"https://github.com/magmasource/alphaMELTS/releases/download/v{version}"

#             # macOS
#             if system == "Darwin":
#                 if "arm" in machine:  # Apple Silicon
#                     url = f"{base_url}/alphamelts-py-{version}-macos-arm64.zip"
#                 else:  # Intel Mac
#                     url = f"{base_url}/alphamelts-py-{version}-macos-x86_64.zip"

#             # Windows
#             elif system == "Windows":
#                 url = f"{base_url}/alphamelts-py-{version}-win64.zip"

#             # Linux (default fallback)
#             elif system == "Linux":
#                 url = f"{base_url}/alphamelts-py-{version}-linux.zip"

#             else:
#                 raise OSError(f"Unsupported system: {system}")
        
#         # if chip == "Apple":
#         #     url = "https://github.com/magmasource/alphaMELTS/releases/download/v2.3.1/alphamelts-py-2.3.1-macos-arm64.zip"
#         # elif chip == "Intel4Mac":
#         #     url = "https://github.com/magmasource/alphaMELTS/releases/download/v2.3.1/alphamelts-py-2.3.1-macos-x86_64.zip"
#         # elif chip == "Windows":
#         #     url = "https://github.com/magmasource/alphaMELTS/releases/download/v2.3.1/alphamelts-py-2.3.1-win64.zip"
#         # else:
#         #     url = "https://github.com/magmasource/alphaMELTS/releases/download/v2.3.1/alphamelts-py-2.3.1-linux.zip"

#         # Path to save the file
#         if file_location is None:
#             if system == "Darwin":
#                 if "arm" in machine:
#                     zip_path = "alphamelts-py-2.3.1-macos-arm64.zip"
#                 else:
#                     zip_path = "alphamelts-py-2.3.1-macos-x86_64.zip"
#             elif system == "Windows":
#                 zip_path = "alphamelts-py-2.3.1-win64.zip"
#             elif system == "Linux":
#                 zip_path = "alphamelts-py-2.3.1-linux.zip"
#             else:
#                 raise OSError(f"Unsupporte system: {system}")
#         else:
#             if system == "Darwin":
#                 if "arm" in machine:
#                     zip_path = file_location + "alphamelts-py-2.3.1-macos-arm64.zip"
#                 else:
#                     zip_path = file_location + "alphamelts-py-2.3.1-macos-x86_64.zip"
#             elif system == "Windows":
#                 zip_path = file_location + "alphamelts-py-2.3.1-win64.zip"
#             elif system == "Linux":
#                 zip_path = file_location + "alphamelts-py-2.3.1-linux.zip"
#             else:
#                 raise OSError(f"Unsupported system: {system}")


#         # Download the file with error handling
#         try:
#             urllib.request.urlretrieve(url, zip_path)
#         except Exception as e:
#             print(f"Error downloading file: {e}. Please contact Matt Gleeson (gleesonm@berkeley.edu) if this error persists.")
#             sys.exit(1)
#             return

#         # Path to extract the contents
#         if file_location is None:
#             extract_path = "alphamelts_py"
#         else:
#             extract_path = file_location + "alphamelts_py"

#         # Extract the zip file with error handling
#         try:
#             with zipfile.ZipFile(zip_path, 'r') as zip_ref:
#                 zip_ref.extractall(extract_path)
#         except zipfile.BadZipFile as e:
#             print(f"Error extracting zip file: {e}")
#             sys.exit(1)
#             return
        
#         if admin is True:
#             # Add the extracted directory to the Python path
#             if file_location is None:
#                 sys.path.append(os.path.join(extract_path, zip_path[:-4]))
#             else:
#                 if system == "Darwin":
#                     if "arm" in machine:
#                         sys.path.append(os.path.join(extract_path,"alphamelts-py-2.3.1-macos-arm64"))
#                     else:
#                         sys.path.append(os.path.join(extract_path,"alphamelts-py-2.3.1-macos-x86_64"))
#                 elif system == "Windows":
#                     sys.path.append(os.path.join(extract_path,"alphamelts-py-2.3.1-win64"))
#                 elif system == "Linux":
#                     sys.path.append(os.path.join(extract_path,"alphamelts-py-2.3.1-linux"))
#                 else:
#                     raise OSError(f"Unsupported system: {system}")

#             try:
#                 import meltsdynamic
#                 importlib.reload(meltsdynamic)
#                 from meltsdynamic import MELTSdynamic
#                 print('Download and Extraction of alphaMELTS for Python files is successful.')
#             except:
#                 print('Error: alphaMELTS for Python not installed correctly.')
#                 return
            
#             # get current working directory
#             cwd = os.getcwd()

#             # Get site-packages path
#             site_packages_path = sysconfig.get_paths()["purelib"]

#             # Define pth file path and custom path
#             pth_file_path = os.path.join(site_packages_path, "my_MELTS_path.pth")
#             if file_location is None:
#                 custom_path = os.path.join(cwd, extract_path, zip_path[:-4])
#             else:
#                 custom_path = os.path.join(zip_path[:-4])

#             # Normalize paths before writing
#             custom_path = os.path.normpath(custom_path)

#             # Write to the .pth file with error handling
#             try:
#                 with open(pth_file_path, "w") as f:
#                     f.write(custom_path)
#             except Exception as e:
#                 print(f"Error writing to .pth file: {e}")
#                 if file_location is None:
#                     custom_path = os.path.join(os.getcwd(), extract_path, zip_path[:-4])
#                 else:
#                     custom_path = os.path.join(zip_path[:-4])
#                 print(f"""The download of the alphaMELTS for Python files was successful, but these were not added to the Python path (likely due to \n absence of administrator privaledges.
#                     The alphaMELTS for Python files are located at 
#                     {custom_path}.
#                     At the start of each notebook using PetThermoTools please add the following lines of code to add these files to the Python path:
#                     import sys
#                     sys.path.append(r"{custom_path}") """)
#                 sys.exit(1)
#         else:
#             if file_location is None:
#                 custom_path = os.path.join(os.getcwd(), extract_path, zip_path[:-4])
#             else:
#                 custom_path = os.path.join(zip_path[:-4])
#             print(f"""alphaMELTS for Python files have been downloaded and are located at 
#                   {custom_path}.
#                   At the start of each notebook using PetThermoTools please add the following lines of code to add these files to the Python path:
#                   import sys
#                   sys.path.append(r"{custom_path}") """)

#         return


def remove_alphaMELTS_path():
    '''
    Removes the custom 'my_MELTS_path.pth' file from the Python site-packages directory.

    This file is what permanently adds the alphaMELTS for Python directory to the Python path.
    Removing it allows you to install or update to a new version without conflicts.

    Note: This function only removes the path file; it does **not** delete the actual
    alphaMELTS files that were previously downloaded. Those must be removed manually.
    '''
    print('Please note this does not remove the previously downloaded alphaMELTS files. That has to be done manually. \n This function simply removes the Python path to those files so that you can install/update to a new version of alphaMELTS for Python.')
    site_packages_dirs = site.getsitepackages()

    for dir in site_packages_dirs:
        print(f"Checking {dir} for .pth files")
        for file in os.listdir(dir):
            if file.endswith(".pth"):
                print(f"Found .pth file: {file}")

    # Path to your site-packages directory (adjust as needed)
    site_packages_path = site.getsitepackages()

    fail = True
    for i in range(len(site_packages_path)):
        # Path to the .pth file
        pth_file_path = os.path.join(site_packages_path[i], "my_MELTS_path.pth")

        # Remove the .pth file
        if os.path.exists(pth_file_path):
            os.remove(pth_file_path)
            fail = False
            print(f"Removed {pth_file_path}")
            break
    
    if fail:
        print(f"Unable to locate {pth_file_path}")

def update_alphaMELTS_path(file_location = None, admin = True, version = "2.3.1"):
    '''
    Updates the installed alphaMELTS for Python files by removing the old path configuration
    and installing a new version.

    This function combines 'remove_alphaMELTS_path' and 'install_alphaMELTS'.

    Args:
        file_location (str, optional): Directory to download and extract the files.
            Passed directly to 'install_alphaMELTS'.
        admin (bool, optional): Controls whether to remove the old path file and
            whether to attempt a permanent path addition for the new version.
            Default is True. Set to False if you don't have admin privileges.
        version (str, optional): The alphaMELTS version to download (e.g., "2.3.2").
            Passed directly to 'install_alphaMELTS'. Default is "2.3.1".

    Returns:
        None. Prints status messages from the underlying removal and installation functions.
    '''
    if admin:
        remove_alphaMELTS_path()

    install_alphaMELTS(file_location = file_location, admin = admin, version = version)

def test_alphaMELTS(path = None):
    '''
    Test the alphaMELTS for Python installation worked! 
    This function should perform a simple fractional crystallization calculation and print the results.
    '''
    if path is not None:
        sys.path.append(path)

    bulk = Compositions['G2']
    
    Res = isobaric_crystallisation(Model = "MELTSv1.0.2",
                                   bulk = bulk,
                                   P_bar = 2000, 
                                   T_end_C = 1100.0,
                                   dt_C = 5.0,
                                   find_liquidus = True,
                                   H2O_init = 0.2)

    print(Res['liquid1'])