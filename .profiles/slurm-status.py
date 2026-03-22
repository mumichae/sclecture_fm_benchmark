#!/usr/bin/env python3
import json
import os
import re
import requests
import subprocess as sp
import shlex
import sys
import time
import logging

logger = logging.getLogger(__name__)

STATUS_ATTEMPTS = 20
SIDECAR_VARS = os.environ.get("SNAKEMAKE_CLUSTER_SIDECAR_VARS", None)
DEBUG = bool(int(os.environ.get("SNAKEMAKE_SLURM_DEBUG", "0")))

if DEBUG:
    logging.basicConfig(level=logging.DEBUG)
    logger.setLevel(logging.DEBUG)


def get_status_direct(jobid):
    """Get status directly from squeue (works without accounting enabled)"""
    for i in range(STATUS_ATTEMPTS):
        try:
            # Use squeue which doesn't require accounting
            squeue_res = sp.check_output(shlex.split(f"squeue -j {jobid} -h -o %T"))
            status = squeue_res.decode().strip()
            if status:
                return status
            else:
                # Job not found in queue (completed or failed)
                return "COMPLETED"
        except sp.CalledProcessError as e:
            logger.debug("squeue returned non-zero: %s", e)
            if i >= STATUS_ATTEMPTS - 1:
                return "COMPLETED"  # Job is done
            else:
                time.sleep(1)
    
    return "COMPLETED"


def get_status_sidecar(jobid):
    """Get status from cluster sidecar"""
    sidecar_vars = json.loads(SIDECAR_VARS)
    url = "http://localhost:%d/job/status/%s" % (sidecar_vars["server_port"], jobid)
    headers = {"Authorization": "Bearer %s" % sidecar_vars["server_secret"]}
    try:
        resp = requests.get(url, headers=headers)
        if resp.status_code == 404:
            return ""  # not found yet
        logger.debug("sidecar returned: %s" % resp.json())
        resp.raise_for_status()
        return resp.json().get("status") or ""
    except requests.exceptions.ConnectionError as e:
        logger.warning("slurm-status.py: could not query side car: %s", e)
        logger.info("slurm-status.py: falling back to direct query")
        return get_status_direct(jobid)


jobid = sys.argv[1]

if SIDECAR_VARS:
    logger.debug("slurm-status.py: querying sidecar")
    status = get_status_sidecar(jobid)
else:
    logger.debug("slurm-status.py: direct query")
    status = get_status_direct(jobid)

logger.debug("job status: %s", repr(status))

if status == "BOOT_FAIL":
    print("failed")
elif status == "OUT_OF_MEMORY":
    print("failed")
elif status.startswith("CANCELLED"):
    print("failed")
elif status == "COMPLETED":
    print("success")
elif status == "DEADLINE":
    print("failed")
elif status == "FAILED":
    print("failed")
elif status == "NODE_FAIL":
    print("failed")
elif status == "PREEMPTED":
    print("failed")
elif status == "TIMEOUT":
    print("failed")
elif status == "SUSPENDED":
    print("running")
else:
    print("running")
