from pydantic import BaseModel
import requests
from xml.etree import ElementTree as ET
from pydantic import BaseModel

class SRR(BaseModel):
    code: str
    num_reads: str | int
    average_reads: dict[int, int]  # index, average reads


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
    
    code: str
    num_reads: int | str
    average_reads: dict[int, int] = {}  # index, average_reads
    for i in tree.iter():
        if i.tag.lower() == "runbundle":
            code = i.attrib["request"]
        if i.tag.lower() == "statistics":
            if "nreads" in i.attrib:  # Collect the number of reads
                if i.attrib["nreads"] == "variable":
                    num_reads = i.attrib["nreads"]
                    break
                else:
                    num_reads = int(i.attrib["nreads"])
        elif i.tag.lower() == "read":
            index: int = int(i.attrib["index"])
            average: int = int(i.attrib["average"])
            average_reads[index] = average
    return SRR(code=code, num_reads=num_reads, average_reads=average_reads)


def collect(srr_code: str) -> SRR:
    url: str = metadata_url(srr_code)
    xml_data = requests_get(url)
    srr: SRR = xml_parser(xml_data)
    return srr


if __name__ == '__main__':
    srr_code: list[str] = [
        "SRR14856518",
        "SRR8387812",
        "SRR20214173"
    ]
    for code in srr_code:
        data = collect(srr_code=code)
        print(data)
