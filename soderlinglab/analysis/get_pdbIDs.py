import re
import time
import json
import zlib
from xml.etree import ElementTree
from urllib.parse import urlparse, parse_qs, urlencode
import requests
from requests.adapters import HTTPAdapter, Retry
import sys
import pandas as pd
import os
import argparse
def check_response(response):
    try:
        response.raise_for_status()
    except requests.HTTPError:
        print(response.json())
        raise


def submit_id_mapping(from_db, to_db, ids):
    request = requests.post(
        f"{API_URL}/idmapping/run",
        data={"from": from_db, "to": to_db, "ids": ",".join(ids),
              }
    )
    check_response(request)
    return request.json()["jobId"]


def get_next_link(headers):
    re_next_link = re.compile(r'<(.+)>; rel="next"')
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)


def check_id_mapping_results_ready(job_id):
    while True:
        request = session.get(f"{API_URL}/idmapping/status/{job_id}")
        check_response(request)
        j = request.json()
        if "jobStatus" in j:
            if j["jobStatus"] == "RUNNING":
                print(f"Retrying in {POLLING_INTERVAL}s")
                time.sleep(POLLING_INTERVAL)
            else:
                raise Exception(j["jobStatus"])
        else:
            return bool(j["results"] or j["failedIds"])


def get_batch(batch_response, file_format, compressed):
    batch_url = get_next_link(batch_response.headers)
    while batch_url:
        batch_response = session.get(batch_url)
        batch_response.raise_for_status()
        yield decode_results(batch_response, file_format, compressed)
        batch_url = get_next_link(batch_response.headers)


def combine_batches(all_results, batch_results, file_format):
    if file_format == "json":
        for key in ("results", "failedIds"):
            if key in batch_results and batch_results[key]:
                all_results[key] += batch_results[key]
    elif file_format == "tsv":
        return all_results + batch_results[1:]
    else:
        return all_results + batch_results
    return all_results


def get_id_mapping_results_link(job_id):
    url = f"{API_URL}/idmapping/details/{job_id}"
    request = session.get(url)
    check_response(request)
    return request.json()["redirectURL"]

def decode_results(response, file_format, compressed):
    if compressed:
        decompressed = zlib.decompress(response.content, 16 + zlib.MAX_WBITS)
        if file_format == "json":
            j = json.loads(decompressed.decode("utf-8"))
            return j
        elif file_format == "tsv":
            return [line for line in decompressed.decode("utf-8").split("\n") if line]
        elif file_format == "xlsx":
            return [decompressed]
        elif file_format == "xml":
            return [decompressed.decode("utf-8")]
        else:
            return decompressed.decode("utf-8")
    elif file_format == "json":
        return response.json()
    elif file_format == "tsv":
        return [line for line in response.text.split("\n") if line]
    elif file_format == "xlsx":
        return [response.content]
    elif file_format == "xml":
        return [response.text]
    return response.text


def get_xml_namespace(element):
    m = re.match(r"\{(.*)\}", element.tag)
    return m.groups()[0] if m else ""


def merge_xml_results(xml_results):
    merged_root = ElementTree.fromstring(xml_results[0])
    for result in xml_results[1:]:
        root = ElementTree.fromstring(result)
        for child in root.findall("{http://uniprot.org/uniprot}entry"):
            merged_root.insert(-1, child)
    ElementTree.register_namespace("", get_xml_namespace(merged_root[0]))
    return ElementTree.tostring(merged_root, encoding="utf-8", xml_declaration=True)


def print_progress_batches(batch_index, size, total):
    n_fetched = min((batch_index + 1) * size, total)
    print(f"Fetched: {n_fetched} / {total}")


def get_id_mapping_results_search(url):
    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    file_format = query["format"][0] if "format" in query else "json"
    if "size" in query:
        size = int(query["size"][0])
    else:
        size = 500
        query["size"] = size
    compressed = (
        query["compressed"][0].lower() == "true" if "compressed" in query else False
    )
    parsed = parsed._replace(query=urlencode(query, doseq=True))
    url = parsed.geturl()
    request = session.get(url)
    check_response(request)
    results = decode_results(request, file_format, compressed)
    total = int(request.headers["x-total-results"])
    print_progress_batches(0, size, total)
    for i, batch in enumerate(get_batch(request, file_format, compressed), 1):
        results = combine_batches(results, batch, file_format)
        print_progress_batches(i, size, total)
    if file_format == "xml":
        return merge_xml_results(results)
    return results


def get_id_mapping_results_stream(url):
    if "/stream/" not in url:
        url = url.replace("/results/", "/results/stream/")
    request = session.get(url)
    check_response(request)
    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    file_format = query["format"][0] if "format" in query else "json"
    compressed = (
        query["compressed"][0].lower() == "true" if "compressed" in query else False
    )
    return decode_results(request, file_format, compressed)

def sort_results(entry):
    check = entry['to']['organism']
    if check.get('taxonId') == 10090:
        # print(entry['to']['primaryAccession'], entry['to']['uniProtkbId'])
        return entry['to']['primaryAccession']

def batch_process_ids(ids, batch_size):
    for i in range(0, len(ids), batch_size):
        yield ids[i:i + batch_size]

def write_results_to_file(results):
    json_path = os.path.expanduser('~/Desktop/D-SCRIPT/data/accession_pdbs.json')
    print(json_path)
    with open(json_path, 'w') as f:
        json.dump(results, f, indent=4)  # Write results in JSON format with indentation for readability
    output_path = os.path.expanduser('~/Desktop/D-SCRIPT/data/pdb_ids.txt')
    with open(output_path, 'w') as f:
        for accid, pdbid in results.items():
            if pdbid is not None:
                 f.write(f'{pdbid},')
def main(GENE_INPUT, taxonID=None):
    global API_URL
    global session
    global POLLING_INTERVAL
    POLLING_INTERVAL = 3
    API_URL = "https://rest.uniprot.org"

    retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))

    GENE_INPUT = [x.split(';')[0] for x in GENE_INPUT]
    batch_size = 1000  # Adjust batch size if needed
    if os.path.isfile(os.path.expanduser('~/Desktop/D-SCRIPT/data/accession_pdbs.json')):
        f = open(os.path.expanduser('~/Desktop/D-SCRIPT/data/accession_pdbs.json'))
        d = json.load(f)
        for acc in GENE_INPUT:
            if acc not in d:
                d[acc] = None
    else:
        d = {g: None for g in GENE_INPUT}

    novelACC = [x for x in d if d[x] is None]
    for batch in batch_process_ids(novelACC, batch_size):
        job_id = submit_id_mapping(from_db="UniProtKB_AC-ID", to_db="PDB", ids=batch)
        if check_id_mapping_results_ready(job_id):
            link = get_id_mapping_results_link(job_id)
            results = get_id_mapping_results_search(link)['results'] #yields list
            for i, entry in enumerate(results):
                if taxonID is not None:
                    d[entry['from']] = sort_results(entry) ## probs need to fix
                else:
                    # print(entry['to'])
                    d[(entry['from'])] = entry["to"]
    write_results_to_file(d)
    return json.dumps(d, indent=4)  # Return dictionary as JSON string

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='get_pdbIDs', description='Map accessions -> pdb')
    parser.add_argument('--ids', type=str, help='Please input the accession pairs file.', required=True)
    parser.add_argument('--taxonID', type=int, help='Is the organism specification needed?')
    parser.add_argument('--sep', type=str, help='How to separate the list of ids?', required=True)
    args = parser.parse_args()
    print(args.sep, file=sys.stderr)
    if args.sep == '\\t':
        args.sep = '\t'
    elif args.sep == '\\n':
        args.sep = '\n'
    if os.path.isfile(args.ids):
        with open(args.ids, 'r') as f:
            GENE_INPUT = [line.strip().split(args.sep) for line in f if line.strip()]
        GENE_INPUT = list(set(item for sublist in GENE_INPUT for item in sublist))
    else:
        GENE_INPUT = args.ids.split(args.sep)

    out = main(GENE_INPUT, args.taxonID)
