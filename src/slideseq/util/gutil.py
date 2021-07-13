# shamelessly heisted from github.com/czbiohub/covidhub

import io
import json
import logging

import pandas as pd
from google.cloud import secretmanager
from google.oauth2.service_account import Credentials
from googleapiclient.discovery import Resource, build

# name for secret storing service account credentials. the API will use existing user credentials
# which should be configured locally to retrieve this secret, and then configure the service account
SECRET_NAME = "projects/velina-208320/secrets/sequencing-credentials/versions/latest"


DIMENSION = "ROWS"
INPUT_MODE = "RAW"
NUM_RETRIES = 5

log = logging.getLogger(__name__)


def get_secrets_manager_credentials(
    secret_name: str = SECRET_NAME,
) -> Credentials:
    client = secretmanager.SecretManagerServiceClient()
    secret_string = client.access_secret_version(name=secret_name).payload.data

    return Credentials.from_service_account_info(json.loads(secret_string))


def get_service(google_credentials: Credentials) -> Resource:
    service = build(
        "drive",
        "v3",
        credentials=google_credentials.with_scopes(
            ["https://www.googleapis.com/auth/drive"]
        ),
        cache_discovery=False,
    )
    return service


class GoogleSheet:
    SHEET_MIMETYPE = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"

    def __init__(self, drive_service: Resource, file_id: str):
        # download the entire spreadsheet as an excel file and store in a buffer
        self.sheet_io = io.BytesIO(
            drive_service.files()
            .export(fileId=file_id, mimeType=GoogleSheet.SHEET_MIMETYPE)
            .execute(num_retries=NUM_RETRIES)
        )

    def __getitem__(self, item):
        """Returns a dataframe from a sheet in the file from its sheet name"""
        return pd.read_excel(self.sheet_io, sheet_name=item)
