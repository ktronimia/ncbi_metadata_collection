#!/usr/bin/env python3

"""
# Dependencies

conda install -c conda-forge ncbi-datasets-cli biopython pandas  requests -y


SELECTION CRITERIA
------------------
1. organism = Organism name
2. genome length > your choice bp
3. EXCLUDE UNVERIFIED records

OUTPUTS
-------
NCBI_collection/
├── NCBI_complete_collection.fasta
├── NCBI_metadata.csv
├── NCBI_accessions.txt
└── GenBank_records/
    └── NCBI_records.gb

"""
print ("\nDear user, Please fillup the following entries\n")

import re
import time
import json
from pathlib import Path
from datetime import datetime

import pandas as pd
from Bio import Entrez
from Bio import SeqIO


# ============================================================
# USER SETTINGS
# ============================================================
print ("\nExample: Human metapneumovirus")
VIRUS_NAME = input("Enter your organism name based on NCBI taxonomy: ").strip()
while True:
    length_input = input("\nMinimum genome length of your organism (in bp): ").strip()
    try:
        MIN_LENGTH = int(length_input)
        break
    except ValueError:
        print("--> Please enter a valid integer for sequence length (e.g., 12000).")

# ------------------------------------------------------------
# IMPORTANT:
# Put your real email here
# ------------------------------------------------------------
print ("\nExample: self@gmail.com")
NCBI_EMAIL = input("Enter your email address (required by NCBI Entrez): ").strip()

# Optional NCBI API key
NCBI_API_KEY = None

# Number of records per NCBI request
BATCH_SIZE = 500

# Output directory
OUTPUT_DIR = Path("NCBI_collection")


# ============================================================
# NCBI SETTINGS
# ============================================================

Entrez.email = NCBI_EMAIL

if NCBI_API_KEY:
    Entrez.api_key = NCBI_API_KEY


# ============================================================
# CREATE DIRECTORIES
# ============================================================

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

GENBANK_DIR = OUTPUT_DIR / "GenBank_records"

GENBANK_DIR.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# OUTPUT FILES
# ============================================================

FASTA_FILE = (
    OUTPUT_DIR /
    "NCBI_complete_collection.fasta"
)

METADATA_FILE = (
    OUTPUT_DIR /
    "NCBI_metadata.csv"
)

ACCESSION_FILE = (
    OUTPUT_DIR /
    "NCBI_accessions.txt"
)

GENBANK_FILE = (
    GENBANK_DIR /
    "NCBI_records.gb"
)


# ============================================================
# SEARCH NCBI
# ============================================================

print("\n" + "=" * 75)
print("NCBI GENOME COLLECTION")
print("=" * 75)

print(f"\nOrganism       : {VIRUS_NAME}")
print(f"Length filter  : > {MIN_LENGTH} bp")
print("Unverified     : EXCLUDED")


# Search organism + length + exclude UNVERIFIED at query level
search_term = (
    f'"{VIRUS_NAME}"[Organism] AND '
    f'{MIN_LENGTH + 1}:1000000[SLEN] NOT '
    f'UNVERIFIED[Title]'
)


print("\nNCBI search:")
print(search_term)


handle = Entrez.esearch(
    db="nuccore",
    term=search_term,
    retmax=0
)

search_result = Entrez.read(handle)

handle.close()


total_records = int(
    search_result["Count"]
)


print(
    f"\nNCBI records matching based on organism + length: "
    f"{total_records}"
)


if total_records == 0:

    raise SystemExit(
        "\nNo records were found. \nPlease check your informations: 'organism name','Taxon ID' 'genome length'"
    )


# ============================================================
# DOWNLOAD ACCESSION NUMBERS
# ============================================================

print("\nDownloading accession list...")


handle = Entrez.esearch(
    db="nuccore",
    term=search_term,
    retstart=0,
    retmax=total_records
)

search_result = Entrez.read(handle)

handle.close()


accessions = list(
    search_result["IdList"]
)


print(
    f"Records retrieved from NCBI: "
    f"{len(accessions)}"
)


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def clean_text(value):
    """
    Clean GenBank qualifier text.
    """

    if value is None:
        return ""

    if isinstance(value, list):

        value = "; ".join(
            str(x)
            for x in value
        )

    value = str(value)

    value = value.replace(
        "\n",
        " "
    )

    value = value.replace(
        "\r",
        " "
    )

    value = re.sub(
        r"\s+",
        " ",
        value
    )

    return value.strip()


# ============================================================
# DATE NORMALIZATION
# ============================================================

def normalize_date(date_string):
    """
    Convert GenBank collection dates to:

        YYYY-MM-DD
    """

    date_string = clean_text(
        date_string
    )

    if not date_string:
        return ""


    date_string = (
        date_string
        .replace("[", "")
        .replace("]", "")
        .strip()
    )


    # --------------------------------------------------------
    # YYYY-MM-DD
    # --------------------------------------------------------

    match = re.fullmatch(
        r"(\d{4})[-/](\d{1,2})[-/](\d{1,2})",
        date_string
    )

    if match:

        year, month, day = match.groups()

        return (
            f"{int(year):04d}-"
            f"{int(month):02d}-"
            f"{int(day):02d}"
        )


    # --------------------------------------------------------
    # YYYY-MM
    # --------------------------------------------------------

    match = re.fullmatch(
        r"(\d{4})[-/](\d{1,2})",
        date_string
    )

    if match:

        year, month = match.groups()

        return (
            f"{int(year):04d}-"
            f"{int(month):02d}-01"
        )


    # --------------------------------------------------------
    # YYYY
    # --------------------------------------------------------

    match = re.fullmatch(
        r"(19\d{2}|20\d{2})",
        date_string
    )

    if match:

        return (
            f"{date_string}-01-01"
        )


    # --------------------------------------------------------
    # Text dates
    # --------------------------------------------------------

    date_formats = [

        "%d-%b-%Y",
        "%d-%B-%Y",

        "%Y-%b-%d",
        "%Y-%B-%d",

        "%d/%m/%Y",
        "%m/%d/%Y",

        "%Y/%m/%d"

    ]


    for fmt in date_formats:

        try:

            date = datetime.strptime(
                date_string,
                fmt
            )

            return date.strftime(
                "%Y-%m-%d"
            )

        except ValueError:

            pass


    # --------------------------------------------------------
    # Extract year
    # --------------------------------------------------------

    year_match = re.search(
        r"(19\d{2}|20\d{2})",
        date_string
    )

    if year_match:

        return (
            f"{year_match.group(1)}-01-01"
        )


    return ""


# ============================================================
# SOURCE FEATURE
# ============================================================

def extract_source_feature(record):

    for feature in record.features:

        if feature.type == "source":

            return feature

    return None


# ============================================================
# UNVERIFIED CHECK
# ============================================================

def is_unverified(record):
    """
    Detect GenBank records marked UNVERIFIED in title, keywords, or comments.
    """
    definition = clean_text(record.description).upper()

    if "UNVERIFIED" in definition:
        return True

    comment = clean_text(record.annotations.get("comment", "")).upper()

    if "UNVERIFIED" in comment:
        return True

    keywords = [k.upper() for k in record.annotations.get("keywords", [])]

    if "UNVERIFIED" in keywords:
        return True

    return False


# ============================================================
# COMPLETE / PARTIAL
# ============================================================

def determine_completeness(record):
    """
    Determine completeness primarily from the sequence title (record.description),
    falling back to keywords and comments if necessary.
    """
    title = clean_text(record.description).lower()

    # Explicit check in title/definition
    if re.search(r"\bcomplete (genome|cds|sequence|gen)\b", title):
        return "complete"
    
    if re.search(r"\bpartial (genome|cds|sequence|gen)\b|\bpartial\b", title):
        return "partial"

    # Fallback to general annotations
    keywords = clean_text(record.annotations.get("keywords", "")).lower()
    comment = clean_text(record.annotations.get("comment", "")).lower()
    full_text = f"{title} {keywords} {comment}"

    if re.search(r"\bcomplete\b", full_text):
        return "complete"
    elif re.search(r"\bpartial\b", full_text):
        return "partial"

    return "unknown"


# ============================================================
# GENOTYPE / SUBTYPE EXTRACTION
# ============================================================

def extract_genotype_and_subtype(
    record,
    source_qualifiers
):

    """
    Extract NCBI genotype/subtype from GenBank metadata.
    """

    genotype = ""
    subtype = ""


    # 1. Direct /genotype qualifier
    if "genotype" in source_qualifiers:

        genotype = clean_text(
            source_qualifiers[
                "genotype"
            ]
        )


    # 2. Direct /subtype qualifier
    if "subtype" in source_qualifiers:

        subtype = clean_text(
            source_qualifiers[
                "subtype"
            ]
        )


    # 3. Search all relevant text
    text_parts = [

        clean_text(
            record.description
        ),

        clean_text(
            record.annotations.get(
                "comment",
                ""
            )
        ),

        clean_text(
            source_qualifiers.get(
                "note",
                ""
            )
        ),

        clean_text(
            source_qualifiers.get(
                "strain",
                ""
            )
        ),

        clean_text(
            source_qualifiers.get(
                "isolate",
                ""
            )
        ),

        clean_text(
            source_qualifiers.get(
                "serotype",
                ""
            )
        )

    ]


    text = " ".join(
        text_parts
    )


    # 4. If genotype explicitly stated
    genotype_patterns = [

        r"\bgenotype\s*[:=]?\s*(A2b1|A2b2|A2a|A2b|A1|A2|B1|B2)\b",

        r"\bgenotype\s+(A2b1|A2b2|A2a|A2b|A1|A2|B1|B2)\b",

        r"\b(NCBI[-_/]?[AB][12](?:a|b)?[12]?)\b"

    ]


    for pattern in genotype_patterns:

        match = re.search(
            pattern,
            text,
            re.IGNORECASE
        )

        if match:

            value = match.group(
                1
            )

            value = re.sub(
                r"^NCBI[-_/]?",
                "",
                value,
                flags=re.IGNORECASE
            )

            if not genotype:

                genotype = value.upper()

            break


    # 5. Detect subtype
    subtype_patterns = [

        r"\bsubtype\s*[:=]?\s*(A2b1|A2b2|A2a|A2b|A1|A2|B1|B2)\b",

        r"\bsub[- ]?genotype\s*[:=]?\s*(A2b1|A2b2|A2a|A2b|A1|A2|B1|B2)\b",

        r"\bsublineage\s*[:=]?\s*(A2b1|A2b2|A2a|A2b|A1|A2|B1|B2)\b"

    ]


    for pattern in subtype_patterns:

        match = re.search(
            pattern,
            text,
            re.IGNORECASE
        )

        if match:

            subtype = (
                match.group(1)
                .upper()
            )

            break


    # 6. Derive broad genotype from subtype if not explicit
    if subtype:

        if subtype.startswith("A"):

            genotype = "A"

        elif subtype.startswith("B"):

            genotype = "B"


    # 7. If only A1/A2/B1/B2 is available
    if not subtype and genotype:

        if genotype.upper() in [
            "A1",
            "A2",
            "B1",
            "B2"
        ]:

            subtype = genotype


    return genotype, subtype


# ============================================================
# NORMALIZED GENBANK METADATA HELPERS
# ============================================================

def qualifier_to_string(value):
    """
    Convert any GenBank qualifier value to one normalized CSV string.
    Multiple values are joined with ' | '.
    """
    if value is None:
        return ""

    if isinstance(value, (list, tuple)):
        cleaned = []
        for item in value:
            item = clean_text(item)
            if item:
                cleaned.append(item)
        return " | ".join(cleaned)

    return clean_text(value)


def get_all_source_qualifiers(records):
    """
    Find every unique qualifier appearing in /source features.

    This makes the CSV schema organism-independent: if one organism
    has /isolation_source, another has /environmental_sample, etc.,
    every observed qualifier becomes its own CSV column.
    """
    qualifier_names = set()

    for record in records:
        source = extract_source_feature(record)
        if source is None:
            continue

        qualifier_names.update(source.qualifiers.keys())

    return sorted(qualifier_names)


def get_all_record_annotation_keys(records):
    """
    Find every GenBank record-level annotation key observed in records.
    """
    keys = set()

    for record in records:
        keys.update(record.annotations.keys())

    return sorted(keys)


def get_feature_summary(record):
    """
    Compact summary of all FEATURES other than /source.

    The complete feature table is still preserved in the .gb file.
    This CSV field is intended as a searchable summary, not a
    replacement for the GenBank FEATURES section.
    """
    feature_parts = []

    for feature in record.features:
        if feature.type == "source":
            continue

        qualifiers = []
        for key, value in feature.qualifiers.items():
            value_string = qualifier_to_string(value)
            if value_string:
                qualifiers.append(f"{key}={value_string}")

        location = clean_text(str(feature.location))

        if qualifiers:
            feature_parts.append(
                f"{feature.type}[{location}]{{" +
                "; ".join(qualifiers) +
                "}"
            )
        else:
            feature_parts.append(
                f"{feature.type}[{location}]"
            )

    return " || ".join(feature_parts)


def get_record_identity_fields(record):
    """
    Extract common GenBank LOCUS/record-level information.
    """
    annotations = record.annotations

    accession = record.id
    base_accession = accession.split(".")[0]

    return {
        "GenBank_accession": accession,
        "GenBank_accession_base": base_accession,
        "GenBank_name": clean_text(record.name),
        "GenBank_definition": clean_text(record.description),
        "GenBank_length": len(record.seq),
        "GenBank_molecule_type": clean_text(
            annotations.get("molecule_type", "")
        ),
        "GenBank_topology": clean_text(
            annotations.get("topology", "")
        ),
        "GenBank_data_file_division": clean_text(
            annotations.get("data_file_division", "")
        ),
        "GenBank_date": clean_text(
            annotations.get("date", "")
        ),
        "GenBank_keywords": qualifier_to_string(
            annotations.get("keywords", "")
        ),
        "GenBank_accessions": qualifier_to_string(
            annotations.get("accessions", "")
        ),
        "GenBank_sequence_version": clean_text(
            annotations.get("sequence_version", "")
        ),
        "GenBank_gi": clean_text(
            annotations.get("gi", "")
        ),
        "GenBank_taxonomy": " | ".join(
            clean_text(x)
            for x in annotations.get("taxonomy", [])
            if clean_text(x)
        ),
        "GenBank_source": clean_text(
            annotations.get("source", "")
        ),
        "GenBank_comment": clean_text(
            annotations.get("comment", "")
        ),
        "GenBank_primary": clean_text(
            annotations.get("primary", "")
        ),
        "GenBank_feature_count": len(record.features),
    }


# ============================================================
# DOWNLOAD GENBANK RECORDS
# ============================================================

print("\n" + "=" * 75)
print("DOWNLOADING GENBANK RECORDS")
print("=" * 75)


all_records = []

verified_records = []

rejected_records = []


for start in range(
    0,
    len(accessions),
    BATCH_SIZE
):

    batch = accessions[
        start:start + BATCH_SIZE
    ]


    print(
        f"\nDownloading records "
        f"{start + 1} - "
        f"{start + len(batch)} "
        f"of {len(accessions)}"
    )


    handle = Entrez.efetch(

        db="nuccore",

        id=batch,

        rettype="gb",

        retmode="text"

    )


    batch_records = list(
        SeqIO.parse(
            handle,
            "genbank"
        )
    )


    handle.close()


    all_records.extend(
        batch_records
    )


    time.sleep(
        0.35
    )


print(
    f"\nGenBank records downloaded: "
    f"{len(all_records)}"
)


# ============================================================
# FILTER RECORDS
# ============================================================

print("\n" + "=" * 75)
print("APPLYING FINAL FILTERS")
print("=" * 75)


for record in all_records:

    # 1. EXCLUDE UNVERIFIED
    if is_unverified(record):

        rejected_records.append(
            (
                record.id,
                "UNVERIFIED"
            )
        )

        continue


    # 2. GENOME LENGTH
    genome_length = len(
        record.seq
    )


    if genome_length <= MIN_LENGTH:

        rejected_records.append(
            (
                record.id,
                "Genome length <= 12000 bp"
            )
        )

        continue


    # 3. SOURCE FEATURE
    source = extract_source_feature(
        record
    )


    if source is None:

        rejected_records.append(
            (
                record.id,
                "No /source feature"
            )
        )

        continue


    qualifiers = source.qualifiers


    # 4. ORGANISM
    organism = clean_text(
        qualifiers.get(
            "organism",
            ""
        )
    )


    if organism.lower() != VIRUS_NAME.lower():

        rejected_records.append(
            (
                record.id,
                f"Organism is '{organism}'"
            )
        )

        continue


    # RECORD PASSED
    verified_records.append(
        record
    )


print(
    f"\nRecords passing ALL filters: "
    f"{len(verified_records)}"
)


print(
    f"Records rejected: "
    f"{len(rejected_records)}"
)


if not verified_records:

    raise SystemExit(
        "\nNo records passed all selection criteria."
    )


# ============================================================
# REMOVE DUPLICATE ACCESSIONS
# ============================================================

unique_records = {}


for record in verified_records:

    accession = record.id

    base_accession = accession.split(
        "."
    )[0]


    if base_accession not in unique_records:

        unique_records[
            base_accession
        ] = record


verified_records = list(
    unique_records.values()
)


print(
    f"Unique records retained: "
    f"{len(verified_records)}"
)


# ============================================================
# SORT BY COLLECTION DATE
# ============================================================

def get_collection_date(record):

    source = extract_source_feature(
        record
    )

    if source is None:
        return ""

    return normalize_date(
        source.qualifiers.get(
            "collection_date",
            ""
        )
    )


verified_records.sort(
    key=lambda r: (
        get_collection_date(r)
        if get_collection_date(r)
        else "9999-99-99"
    )
)


# ============================================================
# WRITE GENBANK
# ============================================================

print(
    "\nWriting GenBank records..."
)


with open(
    GENBANK_FILE,
    "w"
) as handle:

    SeqIO.write(
        verified_records,
        handle,
        "genbank"
    )


# ============================================================
# WRITE FASTA
# ============================================================

print(
    "Writing FASTA..."
)


with open(
    FASTA_FILE,
    "w"
) as handle:

    SeqIO.write(
        verified_records,
        handle,
        "fasta"
    )


# ============================================================
# CREATE NORMALIZED METADATA TABLE
# ============================================================

print(
    "Creating normalized metadata CSV..."
)

# ------------------------------------------------------------
# IMPORTANT:
# Every qualifier observed in ANY /source feature becomes its
# own CSV column. This avoids putting source information into
# one JSON column and keeps the metadata suitable for BEAST,
# R, Python, Excel, etc.
# ------------------------------------------------------------

all_source_qualifiers = get_all_source_qualifiers(
    verified_records
)

all_record_annotations = get_all_record_annotation_keys(
    verified_records
)

print(
    f"Unique /source qualifiers found: "
    f"{len(all_source_qualifiers)}"
)

print(
    "Source qualifiers:",
    ", ".join(all_source_qualifiers)
)

metadata_rows = []

for record in verified_records:

    source = extract_source_feature(record)

    if source is None:
        source_qualifiers = {}
    else:
        source_qualifiers = source.qualifiers

    identity = get_record_identity_fields(record)

    # --------------------------------------------------------
    # Standardized fields
    # --------------------------------------------------------

    organism = qualifier_to_string(
        source_qualifiers.get("organism", "")
    )

    mol_type = qualifier_to_string(
        source_qualifiers.get("mol_type", "")
    )

    isolate = qualifier_to_string(
        source_qualifiers.get("isolate", "")
    )

    strain = qualifier_to_string(
        source_qualifiers.get("strain", "")
    )

    host = qualifier_to_string(
        source_qualifiers.get("host", "")
    )

    geo_loc_name = qualifier_to_string(
        source_qualifiers.get("geo_loc_name", "")
    )

    country = qualifier_to_string(
        source_qualifiers.get("country", "")
    )

    # Keep geo_loc_name as the primary geographic field,
    # but fall back to country when geo_loc_name is absent.
    if not geo_loc_name:
        geo_loc_name = country

    collection_date_raw = qualifier_to_string(
        source_qualifiers.get("collection_date", "")
    )

    collection_date = normalize_date(
        collection_date_raw
    )

    genotype, subtype = extract_genotype_and_subtype(
        record,
        source_qualifiers
    )

    # Preserve the original isolate/strain convenience field.
    if isolate and strain:
        isolate_strain = (
            f"isolate: {isolate}; "
            f"strain: {strain}"
        )
    elif isolate:
        isolate_strain = isolate
    elif strain:
        isolate_strain = strain
    else:
        isolate_strain = ""

    # --------------------------------------------------------
    # Start row with standardized fields
    # --------------------------------------------------------

    row = {
        "accession number": record.id,
        "Genome length": len(record.seq),
        "complete or partial": determine_completeness(record),

        "organism": organism,
        "mol_type": mol_type,
        "isolate or strain": isolate_strain,
        "isolate": isolate,
        "strain": strain,
        "host": host,
        "geo_loc_name": geo_loc_name,
        "country": country,

        "collection_date": collection_date,
        "collection_date_raw": collection_date_raw,

        "genotype": genotype,
        "subtype": subtype,

        # Complete GenBank record-level information.
        **identity,

        # All non-source FEATURES in a compact searchable field.
        "GenBank_feature_summary": get_feature_summary(record),
    }

    # --------------------------------------------------------
    # Add EVERY /source qualifier as its own column.
    #
    # Example:
    # /organism="Human metapneumovirus"
    #       -> source_organism
    #
    # /isolation_source="nasopharyngeal swab"
    #       -> source_isolation_source
    #
    # /collection_date="2019-04-12"
    #       -> source_collection_date
    #
    # This is the key normalization step.
    # --------------------------------------------------------

    for qualifier_name in all_source_qualifiers:

        value = source_qualifiers.get(
            qualifier_name,
            ""
        )

        row[
            f"source_{qualifier_name}"
        ] = qualifier_to_string(value)

    # --------------------------------------------------------
    # Add all record-level annotation keys dynamically.
    # Prefix with GenBank_annotation_ to avoid collisions.
    # --------------------------------------------------------

    for annotation_key in all_record_annotations:

        value = record.annotations.get(
            annotation_key,
            ""
        )

        row[
            f"GenBank_annotation_{annotation_key}"
        ] = qualifier_to_string(value)

    metadata_rows.append(row)


# ------------------------------------------------------------
# Create DataFrame
# ------------------------------------------------------------

metadata = pd.DataFrame(
    metadata_rows
)

# ------------------------------------------------------------
# Normalize missing values to empty strings.
# ------------------------------------------------------------

metadata = metadata.fillna("")


# ------------------------------------------------------------
# Put the most useful BEAST columns first.
# All dynamically discovered source/GenBank columns follow.
# ------------------------------------------------------------

preferred_columns = [
    "accession number",
    "Genome length",
    "complete or partial",
    "organism",
    "mol_type",
    "isolate or strain",
    "isolate",
    "strain",
    "host",
    "geo_loc_name",
    "country",
    "collection_date",
    "collection_date_raw",
    "genotype",
    "subtype",
    "GenBank_definition",
    "GenBank_comment",
    "GenBank_taxonomy",
    "GenBank_keywords",
]

existing_preferred = [
    column
    for column in preferred_columns
    if column in metadata.columns
]

remaining_columns = [
    column
    for column in metadata.columns
    if column not in existing_preferred
]

metadata = metadata[
    existing_preferred +
    remaining_columns
]


# ------------------------------------------------------------
# Save normalized CSV
# ------------------------------------------------------------

metadata.to_csv(
    METADATA_FILE,
    index=False,
    encoding="utf-8-sig"
)

print(
    f"\nNormalized metadata columns: "
    f"{len(metadata.columns)}"
)

print(
    f"Metadata CSV written to:\n  "
    f"{METADATA_FILE}"
)

# ============================================================
# ACCESSION LIST
# ============================================================

with open(
    ACCESSION_FILE,
    "w"
) as handle:

    for accession in metadata[
        "accession number"
    ]:

        handle.write(
            accession + "\n"
        )


# ============================================================
# SUMMARY
# ============================================================

print("\n" + "=" * 75)
print("COLLECTION COMPLETED")
print("=" * 75)

print(f"\nVirus                 : {VIRUS_NAME}")
print(f"Length requirement     : > {MIN_LENGTH} bp")
print(f"Unverified sequences   : EXCLUDED")
print(f"\nRecords found          : {total_records}")
print(f"Records downloaded     : {len(all_records)}")
print(f"Records rejected       : {len(rejected_records)}")
print(f"Records retained       : {len(metadata)}")
print(f"\nFASTA:\n  {FASTA_FILE}")
print(f"\nMetadata CSV:\n  {METADATA_FILE}")
print(f"\nAccession list:\n  {ACCESSION_FILE}")
print(f"\nGenBank:\n  {GENBANK_FILE}")

print("\nFirst few records:\n")
print(metadata.head(10).to_string(index=False))
print("\nDone.")
