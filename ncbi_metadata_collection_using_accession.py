#!/usr/bin/env python3

"""

# Dependencies

conda install -c conda-forge ncbi-datasets-cli biopython pandas  requests -y

# Create "accession_list.txt" in the same directory in the following format
accession1
accession2
accession3
...

    
"""


# ============================================================
# IMPORTS
# ============================================================

import re
import time
from pathlib import Path
from datetime import datetime

import pandas as pd
from Bio import Entrez
from Bio import SeqIO


# ============================================================
# USER SETTINGS
# ============================================================

# ------------------------------------------------------------

ACCESSION_INPUT_FILE = Path("accession_list.txt")


# ------------------------------------------------------------
# OUTPUT DIRECTORY
# ------------------------------------------------------------

OUTPUT_DIR = Path("GenBank_collection")


# ------------------------------------------------------------
# NCBI EMAIL
#
# Required by NCBI Entrez.
# Replace with your own email.
# ------------------------------------------------------------

NCBI_EMAIL = "your_email@example.com"


# ------------------------------------------------------------
# OPTIONAL NCBI API KEY
#
# Keep None if you do not have one.
# ------------------------------------------------------------

NCBI_API_KEY = None


# ------------------------------------------------------------
# DOWNLOAD BATCH SIZE
# ------------------------------------------------------------

BATCH_SIZE = 500

MIN_LENGTH = None
ORGANISM_FILTER = None
TAXON_ID = None

EXCLUDE_UNVERIFIED = True
REQUIRE_SOURCE = True


# ------------------------------------------------------------
# NORMALIZATION SETTINGS
# ------------------------------------------------------------

NORMALIZE_COLLECTION_DATE = True


# ============================================================
# NCBI CONFIGURATION
# ============================================================

Entrez.email = NCBI_EMAIL

if NCBI_API_KEY:
    Entrez.api_key = NCBI_API_KEY


# ============================================================
# OUTPUT FILES
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


FASTA_FILE = OUTPUT_DIR / "complete_collection.fasta"

METADATA_FILE = OUTPUT_DIR / "metadata.csv"

ACCESSION_FILE = OUTPUT_DIR / "accessions.txt"

GENBANK_FILE = GENBANK_DIR / "records.gb"


# ============================================================
# BASIC HELPERS
# ============================================================

def clean_text(value):
    """
    Convert GenBank values to clean, single-line text.
    """

    if value is None:
        return ""

    if isinstance(value, (list, tuple)):
        value = " | ".join(
            str(x)
            for x in value
            if x is not None
        )

    value = str(value)

    value = value.replace("\n", " ")
    value = value.replace("\r", " ")
    value = re.sub(r"\s+", " ", value)

    return value.strip()


def qualifier_to_string(value):
    """
    Convert a GenBank qualifier into a normalized CSV value.

    Multiple values are retained and separated by ' | '.
    """

    if value is None:
        return ""

    if isinstance(value, (list, tuple)):
        values = []

        for item in value:
            item = clean_text(item)

            if item:
                values.append(item)

        return " | ".join(values)

    return clean_text(value)


def safe_column_name(value):
    """
    Convert arbitrary GenBank qualifier names into safe CSV column names.
    """

    value = clean_text(value)

    value = re.sub(
        r"[^A-Za-z0-9_.-]+",
        "_",
        value
    )

    value = value.strip("_")

    return value or "unknown"


# ============================================================
# DATE NORMALIZATION
# ============================================================

def normalize_date(date_string):
    """
    Normalize common GenBank dates to:

        YYYY-MM-DD

    Year-only dates become YYYY-01-01.
    Year-month dates become YYYY-MM-01.

    If the value cannot be interpreted, the original cleaned
    value is returned rather than silently deleting information.
    """

    date_string = clean_text(date_string)

    if not date_string:
        return ""

    date_string = (
        date_string
        .replace("[", "")
        .replace("]", "")
        .strip()
    )

    # YYYY-MM-DD
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

    # YYYY-MM
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

    # YYYY
    match = re.fullmatch(
        r"(19\d{2}|20\d{2})",
        date_string
    )

    if match:
        return f"{date_string}-01-01"

    # Text dates
    date_formats = [
        "%d-%b-%Y",
        "%d-%B-%Y",
        "%Y-%b-%d",
        "%Y-%B-%d",
        "%d/%m/%Y",
        "%m/%d/%Y",
        "%Y/%m/%d",
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

    # If there is a year, preserve the information.
    year_match = re.search(
        r"(19\d{2}|20\d{2})",
        date_string
    )

    if year_match:
        return f"{year_match.group(1)}-01-01"

    # Never destroy an unrecognized date.
    return date_string


# ============================================================
# ACCESSION INPUT
# ============================================================

def read_accession_list(path):
    """
    Read accession numbers from a plain-text accession list.
    """

    if not path.exists():

        raise FileNotFoundError(
            f"Accession list not found:\n{path}"
        )

    accessions = []

    with open(
        path,
        "r",
        encoding="utf-8"
    ) as handle:

        for line in handle:

            line = line.strip()

            if not line:
                continue

            if line.startswith("#"):
                continue

            # Allow whitespace-separated accession lists too.
            parts = line.split()

            for accession in parts:

                accession = accession.strip()

                if accession:
                    accessions.append(
                        accession
                    )

    # Remove duplicates while preserving order.
    unique = []

    seen = set()

    for accession in accessions:

        base = accession.split(".")[0]

        if base not in seen:

            seen.add(base)

            unique.append(
                accession
            )

    return unique


# ============================================================
# SOURCE FEATURE
# ============================================================

def extract_source_features(record):
    """
    Return all /source features.

    Normally there is one source feature per record, but the function
    intentionally supports multiple source features.
    """

    return [
        feature
        for feature in record.features
        if feature.type == "source"
    ]


# ============================================================
# UNVERIFIED CHECK
# ============================================================

def is_unverified(record):

    definition = clean_text(
        record.description
    ).upper()

    if "UNVERIFIED" in definition:
        return True

    comment = clean_text(
        record.annotations.get(
            "comment",
            ""
        )
    ).upper()

    if "UNVERIFIED" in comment:
        return True

    keywords = record.annotations.get(
        "keywords",
        []
    )

    if isinstance(keywords, str):
        keywords = [keywords]

    for keyword in keywords:

        if "UNVERIFIED" in clean_text(
            keyword
        ).upper():

            return True

    return False


# ============================================================
# COMPLETENESS
# ============================================================

def determine_completeness(record):

    text = " ".join([
        clean_text(record.description),
        clean_text(
            record.annotations.get(
                "comment",
                ""
            )
        ),
        clean_text(
            record.annotations.get(
                "keywords",
                ""
            )
        ),
    ]).lower()

    if re.search(
        r"\bcomplete\s+(genome|sequence|cds|segment)\b",
        text
    ):
        return "complete"

    if re.search(
        r"\bpartial\s+(genome|sequence|cds|segment)\b",
        text
    ):
        return "partial"

    if re.search(
        r"\bpartial\b",
        text
    ):
        return "partial"

    if re.search(
        r"\bcomplete\b",
        text
    ):
        return "complete"

    return "unknown"


# ============================================================
# RECORD-LEVEL GENBANK INFORMATION
# ============================================================

def get_record_annotations(record):
    """
    Dynamically extract all record-level GenBank annotations.

    Nothing is restricted to a predefined organism.
    """

    data = {}

    for key, value in record.annotations.items():

        column = (
            "GenBank_annotation_"
            + safe_column_name(key)
        )

        data[column] = qualifier_to_string(
            value
        )

    return data


def get_record_identity(record):
    """
    Standard GenBank record information.
    """

    accession = clean_text(
        record.id
    )

    base_accession = accession.split(".")[0]

    return {
        "GenBank_accession":
            accession,

        "GenBank_accession_base":
            base_accession,

        "GenBank_name":
            clean_text(record.name),

        "GenBank_definition":
            clean_text(record.description),

        "GenBank_length":
            len(record.seq),

        "GenBank_feature_count":
            len(record.features),
    }


# ============================================================
# ALL /SOURCE INFORMATION
# ============================================================

def collect_source_qualifiers(records):
    """
    Discover EVERY unique qualifier appearing in ANY /source feature.
    """

    qualifiers = set()

    for record in records:

        for source in extract_source_features(
            record
        ):

            qualifiers.update(
                source.qualifiers.keys()
            )

    return sorted(
        qualifiers
    )


def source_values_by_qualifier(record):
    """
    Collect values from all /source features.

    If a record has multiple /source features, their values are
    combined with ' | '.

    This prevents loss of source information.
    """

    result = {}

    for source in extract_source_features(
        record
    ):

        for key, value in source.qualifiers.items():

            if key not in result:
                result[key] = []

            values = value

            if not isinstance(
                values,
                (list, tuple)
            ):
                values = [values]

            for item in values:

                item = clean_text(item)

                if item:
                    result[key].append(
                        item
                    )

    return {
        key: " | ".join(
            values
        )
        for key, values in result.items()
    }


def source_raw_text(record):
    """
    Preserve all /source qualifiers as a readable raw-style string.

    This provides an additional safety net alongside the individual
    source_* columns.
    """

    blocks = []

    for source in extract_source_features(
        record
    ):

        parts = []

        for key, value in source.qualifiers.items():

            values = value

            if not isinstance(
                values,
                (list, tuple)
            ):
                values = [values]

            for item in values:

                item = clean_text(item)

                if item:

                    parts.append(
                        f"/{key}={item}"
                    )

        if parts:

            blocks.append(
                "; ".join(parts)
            )

    return " || ".join(
        blocks
    )


# ============================================================
# ALL FEATURE QUALIFIERS
# ============================================================

def collect_feature_columns(records):
    """
    Discover all feature-type/qualifier combinations.

    Example:

        CDS + gene
            -> feature_CDS_gene

        CDS + product
            -> feature_CDS_product

        gene + gene
            -> feature_gene_gene

    Values from repeated features are retained.
    """

    combinations = set()

    for record in records:

        for feature in record.features:

            feature_type = safe_column_name(
                feature.type
            )

            for qualifier in feature.qualifiers:

                qualifier_name = safe_column_name(
                    qualifier
                )

                combinations.add(
                    (
                        feature_type,
                        qualifier_name
                    )
                )

    return sorted(
        combinations
    )


def get_feature_metadata(record):
    """
    Extract all feature qualifiers into normalized aggregate columns.

    The .gb file remains the authoritative representation of exact
    feature locations and structure.
    """

    data = {}

    for feature in record.features:

        feature_type = safe_column_name(
            feature.type
        )

        for qualifier, value in feature.qualifiers.items():

            qualifier_name = safe_column_name(
                qualifier
            )

            column = (
                f"feature_{feature_type}_"
                f"{qualifier_name}"
            )

            value_string = qualifier_to_string(
                value
            )

            if not value_string:
                continue

            if column not in data:
                data[column] = []

            data[column].append(
                value_string
            )

    for column in data:

        # Preserve repeated feature values without duplicates.
        values = []

        for value in data[column]:

            if value not in values:
                values.append(value)

        data[column] = " | ".join(
            values
        )

    return data


def get_feature_summary(record):
    """
    Human-readable summary of every feature.

    Includes feature type, location and qualifiers.
    """

    feature_parts = []

    for feature in record.features:

        feature_type = clean_text(
            feature.type
        )

        location = clean_text(
            str(feature.location)
        )

        qualifiers = []

        for key, value in feature.qualifiers.items():

            value_string = qualifier_to_string(
                value
            )

            if value_string:

                qualifiers.append(
                    f"{key}={value_string}"
                )

        if qualifiers:

            feature_parts.append(
                f"{feature_type}[{location}]{{"
                + "; ".join(qualifiers)
                + "}"
            )

        else:

            feature_parts.append(
                f"{feature_type}[{location}]"
            )

    return " || ".join(
        feature_parts
    )


# ============================================================
# FILTERING
# ============================================================

def record_passes_filters(record):

    # --------------------------------------------------------
    # UNVERIFIED
    # --------------------------------------------------------

    if (
        EXCLUDE_UNVERIFIED
        and is_unverified(record)
    ):
        return False, "UNVERIFIED"


    # --------------------------------------------------------
    # LENGTH
    # --------------------------------------------------------

    if (
        MIN_LENGTH is not None
        and len(record.seq) < MIN_LENGTH
    ):

        return (
            False,
            f"length < {MIN_LENGTH}"
        )


    # --------------------------------------------------------
    # ORGANISM
    # --------------------------------------------------------

    if ORGANISM_FILTER:

        source_features = (
            extract_source_features(record)
        )

        organisms = []

        for source in source_features:

            organisms.extend(
                source.qualifiers.get(
                    "organism",
                    []
                )
                if isinstance(
                    source.qualifiers.get(
                        "organism",
                        []
                    ),
                    list
                )
                else [
                    source.qualifiers.get(
                        "organism",
                        ""
                    )
                ]
            )

        organisms = [
            clean_text(x).lower()
            for x in organisms
            if clean_text(x)
        ]

        if (
            ORGANISM_FILTER.lower()
            not in organisms
        ):

            return (
                False,
                "organism filter mismatch"
            )


    # --------------------------------------------------------
    # TAXON ID
    # --------------------------------------------------------

    if TAXON_ID:

        taxon_target = (
            f"taxon:{TAXON_ID}".lower()
        )

        found = False

        for source in extract_source_features(
            record
        ):

            values = source.qualifiers.get(
                "db_xref",
                []
            )

            if not isinstance(
                values,
                (list, tuple)
            ):
                values = [values]

            for value in values:

                if (
                    clean_text(value).lower()
                    == taxon_target
                ):

                    found = True
                    break

            if found:
                break

        if not found:

            return (
                False,
                f"taxon:{TAXON_ID} not found"
            )


    # --------------------------------------------------------
    # SOURCE FEATURE
    # --------------------------------------------------------

    if (
        REQUIRE_SOURCE
        and not extract_source_features(
            record
        )
    ):

        return (
            False,
            "no /source feature"
        )


    return True, "passed"


# ============================================================
# DOWNLOAD GENBANK RECORDS
# ============================================================

def download_records(accessions):

    all_records = []

    print("\n" + "=" * 75)
    print("DOWNLOADING GENBANK RECORDS")
    print("=" * 75)

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

        try:

            batch_records = list(
                SeqIO.parse(
                    handle,
                    "genbank"
                )
            )

        finally:

            handle.close()

        all_records.extend(
            batch_records
        )

        time.sleep(0.35)

    return all_records


# ============================================================
# CREATE METADATA CSV
# ============================================================

def create_metadata(records):

    print(
        "\nCreating normalized metadata CSV..."
    )

    # Discover schema dynamically from all retained records.
    source_qualifiers = (
        collect_source_qualifiers(records)
    )

    feature_columns = (
        collect_feature_columns(records)
    )

    print(
        f"Unique /source qualifiers: "
        f"{len(source_qualifiers)}"
    )

    if source_qualifiers:

        print(
            "Source qualifiers:"
        )

        print(
            "  "
            + ", ".join(
                source_qualifiers
            )
        )

    print(
        f"Feature type/qualifier columns: "
        f"{len(feature_columns)}"
    )


    rows = []


    for record in records:

        source_values = (
            source_values_by_qualifier(
                record
            )
        )

        identity = (
            get_record_identity(
                record
            )
        )

        row = {

            # ------------------------------------------------
            # Standardized general fields
            # ------------------------------------------------

            "accession number":
                record.id,

            "Genome length":
                len(record.seq),

            "complete or partial":
                determine_completeness(
                    record
                ),

            "organism":
                source_values.get(
                    "organism",
                    ""
                ),

            "mol_type":
                source_values.get(
                    "mol_type",
                    ""
                ),

            "isolate":
                source_values.get(
                    "isolate",
                    ""
                ),

            "strain":
                source_values.get(
                    "strain",
                    ""
                ),

            "isolate or strain":
                " | ".join(
                    value
                    for value in [
                        source_values.get(
                            "isolate",
                            ""
                        ),
                        source_values.get(
                            "strain",
                            ""
                        ),
                    ]
                    if value
                ),

            "host":
                source_values.get(
                    "host",
                    ""
                ),

            "geo_loc_name":
                source_values.get(
                    "geo_loc_name",
                    ""
                ),

            "country":
                source_values.get(
                    "country",
                    ""
                ),

            "collection_date_raw":
                source_values.get(
                    "collection_date",
                    ""
                ),

            "collection_date":
                (
                    normalize_date(
                        source_values.get(
                            "collection_date",
                            ""
                        )
                    )
                    if NORMALIZE_COLLECTION_DATE
                    else source_values.get(
                        "collection_date",
                        ""
                    )
                ),

            # ------------------------------------------------
            # Generic genotype/subtype fields.
            #
            # NO organism-specific interpretation is applied.
            # If /genotype or /subtype exists, it is copied.
            # ------------------------------------------------

            "genotype":
                source_values.get(
                    "genotype",
                    ""
                ),

            "subtype":
                source_values.get(
                    "subtype",
                    ""
                ),

            # ------------------------------------------------
            # GenBank identity
            # ------------------------------------------------

            **identity,

            # ------------------------------------------------
            # All /source information as one backup field
            # ------------------------------------------------

            "GenBank_source_qualifiers":
                source_raw_text(
                    record
                ),

            # ------------------------------------------------
            # Complete feature summary
            # ------------------------------------------------

            "GenBank_feature_summary":
                get_feature_summary(
                    record
                ),
        }


        # ----------------------------------------------------
        # EVERY /source QUALIFIER -> source_<qualifier>
        # ----------------------------------------------------

        for qualifier in source_qualifiers:

            column = (
                "source_"
                + safe_column_name(
                    qualifier
                )
            )

            row[column] = source_values.get(
                qualifier,
                ""
            )


        # ----------------------------------------------------
        # ALL RECORD-LEVEL ANNOTATIONS
        # ----------------------------------------------------

        row.update(
            get_record_annotations(
                record
            )
        )


        # ----------------------------------------------------
        # ALL FEATURE QUALIFIERS
        # ----------------------------------------------------

        feature_data = (
            get_feature_metadata(
                record
            )
        )

        for feature_type, qualifier in feature_columns:

            column = (
                f"feature_{feature_type}_"
                f"{qualifier}"
            )

            row[column] = feature_data.get(
                column,
                ""
            )


        rows.append(row)


    metadata = pd.DataFrame(
        rows
    )

    metadata = metadata.fillna("")


    # --------------------------------------------------------
    # Preferred columns first.
    # Dynamic columns follow.
    # --------------------------------------------------------

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

        "GenBank_accession",
        "GenBank_accession_base",
        "GenBank_name",
        "GenBank_definition",
        "GenBank_length",
        "GenBank_feature_count",

        "GenBank_source_qualifiers",
        "GenBank_feature_summary",
    ]


    existing = [
        column
        for column in preferred_columns
        if column in metadata.columns
    ]


    remaining = [
        column
        for column in metadata.columns
        if column not in existing
    ]


    metadata = metadata[
        existing + remaining
    ]


    metadata.to_csv(
        METADATA_FILE,
        index=False,
        encoding="utf-8-sig"
    )


    print(
        f"\nMetadata columns created: "
        f"{len(metadata.columns)}"
    )

    print(
        f"Metadata CSV:\n  "
        f"{METADATA_FILE}"
    )

    return metadata


# ============================================================
# MAIN
# ============================================================

def main():

    print("\n" + "=" * 75)
    print("UNIVERSAL NCBI GENBANK COLLECTION")
    print("=" * 75)

    print(
        "\nInput accession list:"
        f" {ACCESSION_INPUT_FILE}"
    )

    print(
        f"Output directory:"
        f" {OUTPUT_DIR}"
    )

    print(
        f"Minimum length:"
        f" {MIN_LENGTH if MIN_LENGTH is not None else 'NONE'}"
    )

    print(
        f"Organism filter:"
        f" {ORGANISM_FILTER if ORGANISM_FILTER else 'NONE'}"
    )

    print(
        f"Taxon filter:"
        f" {TAXON_ID if TAXON_ID else 'NONE'}"
    )

    print(
        f"Exclude UNVERIFIED:"
        f" {EXCLUDE_UNVERIFIED}"
    )

    print(
        f"Require /source:"
        f" {REQUIRE_SOURCE}"
    )


    # --------------------------------------------------------
    # Read accession list
    # --------------------------------------------------------

    accessions = read_accession_list(
        ACCESSION_INPUT_FILE
    )

    if not accessions:

        raise SystemExit(
            "\nNo accession numbers were found."
        )

    print(
        f"\nUnique accessions:"
        f" {len(accessions)}"
    )


    # --------------------------------------------------------
    # Download
    # --------------------------------------------------------

    all_records = download_records(
        accessions
    )

    print(
        f"\nGenBank records downloaded:"
        f" {len(all_records)}"
    )


    if not all_records:

        raise SystemExit(
            "\nNo GenBank records were downloaded."
        )


    # --------------------------------------------------------
    # Apply filters
    # --------------------------------------------------------

    print("\n" + "=" * 75)
    print("APPLYING FILTERS")
    print("=" * 75)

    retained = []
    rejected = []

    for record in all_records:

        passed, reason = (
            record_passes_filters(
                record
            )
        )

        if passed:

            retained.append(
                record
            )

        else:

            rejected.append(
                (
                    record.id,
                    reason
                )
            )


    print(
        f"\nRecords retained:"
        f" {len(retained)}"
    )

    print(
        f"Records rejected:"
        f" {len(rejected)}"
    )


    if not retained:

        raise SystemExit(
            "\nNo records passed the selected filters."
        )


    # --------------------------------------------------------
    # Remove duplicate accession versions
    # --------------------------------------------------------

    unique_records = {}

    for record in retained:

        base = record.id.split(".")[0]

        if base not in unique_records:

            unique_records[
                base
            ] = record


    retained = list(
        unique_records.values()
    )


    print(
        f"Unique records retained:"
        f" {len(retained)}"
    )


    # --------------------------------------------------------
    # Sort by collection date when available.
    # Undated records go last.
    # --------------------------------------------------------

    def sort_key(record):

        source_values = (
            source_values_by_qualifier(
                record
            )
        )

        date = source_values.get(
            "collection_date",
            ""
        )

        if NORMALIZE_COLLECTION_DATE:

            date = normalize_date(
                date
            )

        return (
            date
            if date
            else "9999-99-99"
        )


    retained.sort(
        key=sort_key
    )


    # --------------------------------------------------------
    # Write GenBank
    # --------------------------------------------------------

    print(
        "\nWriting GenBank file..."
    )

    with open(
        GENBANK_FILE,
        "w",
        encoding="utf-8"
    ) as handle:

        SeqIO.write(
            retained,
            handle,
            "genbank"
        )


    # --------------------------------------------------------
    # Write FASTA
    # --------------------------------------------------------

    print(
        "Writing FASTA..."
    )

    with open(
        FASTA_FILE,
        "w",
        encoding="utf-8"
    ) as handle:

        SeqIO.write(
            retained,
            handle,
            "fasta"
        )


    # --------------------------------------------------------
    # Create metadata
    # --------------------------------------------------------

    metadata = create_metadata(
        retained
    )


    # --------------------------------------------------------
    # Write final accession list
    # --------------------------------------------------------

    with open(
        ACCESSION_FILE,
        "w",
        encoding="utf-8"
    ) as handle:

        for accession in metadata[
            "accession number"
        ]:

            handle.write(
                accession + "\n"
            )


    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    print("\n" + "=" * 75)
    print("COLLECTION COMPLETED")
    print("=" * 75)

    print(
        f"\nInput accessions       : "
        f"{len(accessions)}"
    )

    print(
        f"Downloaded records     : "
        f"{len(all_records)}"
    )

    print(
        f"Rejected records       : "
        f"{len(rejected)}"
    )

    print(
        f"Retained records       : "
        f"{len(retained)}"
    )

    print(
        f"CSV columns            : "
        f"{len(metadata.columns)}"
    )

    print(
        f"\nFASTA:\n  {FASTA_FILE}"
    )

    print(
        f"\nMetadata CSV:\n  {METADATA_FILE}"
    )

    print(
        f"\nAccession list:\n  {ACCESSION_FILE}"
    )

    print(
        f"\nGenBank:\n  {GENBANK_FILE}"
    )

    if rejected:

        print(
            "\nRejected records:"
        )

        for accession, reason in rejected[:20]:

            print(
                f"  {accession}: {reason}"
            )

        if len(rejected) > 20:

            print(
                f"  ... and "
                f"{len(rejected) - 20} more"
            )

    print("\nDone.")


# ============================================================
# RUN
# ============================================================

if __name__ == "__main__":
    main()
