import subprocess

def download_uniprot(accession: str) -> str:
    """Download record from uniprot using
    uniprot rest api

    Args:
        accession: requesting accession
    Returns:
        str(uniprot record)
    """

    command = f'curl -X GET "https://rest.uniprot.org/uniprotkb/{accession}.fasta" -H "accept: text/plain"'
    uniprot_record = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL).stdout.read().decode()

    return uniprot_record

def read_uniprot(uniprot_record:str) -> dict:
    """read record from uniprot download 
    request

    Args:
        uniprot_record: subject uniprot record
    Returns:
        dict("header":str(), "sequence":str())
    """

    uniprot_info = dict()

    for line in uniprot_record.splitlines():
        if line.startswith(">"):
            uniprot_info["header"] = line.strip().strip(">")
        else:
            try:
                uniprot_info["sequence"] += line.strip().upper()
            except:
                uniprot_info["sequence"] = line.upper()
    
    return uniprot_info

def parse_uniprot(accession:str) -> dict:
    """parse record from given accession

    Args:
        uniprot_record: subject uniprot record
    Returns:
        dict("header":str(), "sequence":str())
    """

    uniprot_info = download_uniprot(accession)
    uniprot_info = read_uniprot(uniprot_info)
    
    return uniprot_info