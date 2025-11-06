import urllib.request
import zipfile
import os
import sys
import site
import sysconfig
import importlib
from pathlib import Path

def install_MAGEMinCalc():
    '''
    Establish a new julia environment that will be used for any MAGEMin calculations performed through PetThermoTools.
    MAGEMinCalc and MAGEMin_C are added to this new environment.
    '''
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
            Pkg.add(url="https://github.com/gleesonm1/MAGEMinCalc.git", rev="v0.5.0")
        catch e
            @warn "Failed to install MAGEMinCalc via HTTPS, retrying..." exception=e
            ENV["JULIA_SSL_CA_ROOTS_PATH"] = ""
            Pkg.add(url="https://github.com/gleesonm1/MAGEMinCalc.git", rev="v0.5.0")
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
            Pkg.add(url="https://github.com/gleesonm1/MAGEMinCalc.git", rev="v0.5.0")
        catch e
            @warn "Failed to install MAGEMinCalc via HTTPS, retrying..." exception=e
            ENV["JULIA_SSL_CA_ROOTS_PATH"] = ""
            Pkg.add(url="https://github.com/gleesonm1/MAGEMinCalc.git", rev="v0.5.0")
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

def install_alphaMELTS(chip="Linux", file_location = None, admin = False):
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
            print(f"Error downloading file: {e}. Please contact Matt Gleeson (gleesonm@berkeley.edu) if this error persists.")
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
        
        if admin is True:
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
                if file_location is None:
                    custom_path = os.path.join(os.getcwd(), extract_path, zip_path[:-4])
                else:
                    custom_path = os.path.join(zip_path[:-4])
                print(f"""The download of the alphaMELTS for Python files was successful, but these were not added to the Python path (likely due to \n absence of administrator privaledges.
                    The alphaMELTS for Python files are located at 
                    {custom_path}.
                    At the start of each notebook using PetThermoTools please add the following lines of code to add these files to the Python path:
                    import sys
                    sys.path.append(r"{custom_path}") """)
                sys.exit(1)
        else:
            if file_location is None:
                custom_path = os.path.join(os.getcwd(), extract_path, zip_path[:-4])
            else:
                custom_path = os.path.join(zip_path[:-4])
            print(f"""alphaMELTS for Python files have been downloaded and are located at 
                  {custom_path}.
                  At the start of each notebook using PetThermoTools please add the following lines of code to add these files to the Python path:
                  import sys
                  sys.path.append(r"{custom_path}") """)

        return
    
def update_alphaMELTS_path(chip="Linux", file_location = None, admin = False):
    '''
    Update the alphaMELTS for Python files to the latest version. \n Please keep an eye of the alphaMELTS Discord Server for information on when a new version of alphaMELTS for Python is published.
    '''
    if admin:
        print('Please note this does not remove the previously downloaded alphaMELTS files. That has to be done manually. \n This function simply removes the Python path to those files so that you can install/update to a new version of alphaMELTS for Python.')
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

    install_alphaMELTS(chip=chip, file_location = file_location, admin = admin)

    return