from pydantic import BaseModel
import requests
from xml.etree import ElementTree as ET
from pydantic import BaseModel

class SRR(BaseModel):
    code: str = ""
    num_reads: int = -1
    average_reads: dict[int, int] = {} # index, average reads
    I_file_index: int = -1
    R1_file_index: int = -1
    R2_file_index: int = -1


def metadata_url(srr_code: str) -> str:
    return f"https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/run_new?acc={srr_code}"


def requests_get(url: str):
    """
    This function will download the webpage using requests, which can tell us various details about the page in the "loading" state
    :param url:
    :return:
    """
    request = requests.get(url)
    return request.text


def xml_parser(xml_data) -> SRR:
    root = ET.fromstring(xml_data)
    tree = ET.ElementTree(element=root)
    
    srr_builder = SRR()
    filename_index: int = 1
    for i in tree.iter():
        #print(i.tag, i.attrib)
        match i.tag.lower():
            case "runbundle":
                srr_builder.code = i.attrib["request"]
            case "statistics":
                if "nreads" in i.attrib:  # Collect the number of reads
                    if i.attrib["nreads"] == "variable":
                        srr_builder.num_reads = 2
                    else:
                        srr_builder.num_reads = i.attrib["nreads"]
            case "read":
                index: int = int(i.attrib["index"])
                average: int = int(i.attrib["average"])
                srr_builder.average_reads[index] = average
            case "srafile":
                filename = i.attrib["filename"]
                if "I1" in filename:
                    srr_builder.I_file_index = filename_index
                elif "R1" in filename:
                    srr_builder.R1_file_index = filename_index
                elif "R2" in filename and filename[0:4] != "SRR2":  # Make sure we dont get an SRR code like "SRR23456"
                    srr_builder.R2_file_index = filename_index
                filename_index += 1
            
    return srr_builder


def collect(srr_code: str) -> SRR:
    url: str = metadata_url(srr_code)
    xml_data = requests_get(url)
    srr: SRR = xml_parser(xml_data)
    return srr


if __name__ == '__main__':
    srr_code: list[str] = [
        # "SRR22474175",  # Two reads
        # "SRR8387812",  # Three reads
        "SRR20214173"   # Variable reads
    ]
    for code in srr_code:
        data = collect(srr_code=code)
        print(data)
