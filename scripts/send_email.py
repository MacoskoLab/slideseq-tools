#!/usr/bin/python

# This script is to send an email to the user(s)

import smtplib
import sys

if len(sys.argv) != 4:
    print("Please provide three arguments: receiver, subject, and email content!")
    sys.exit()

FROM = "slideseq@gmail.com"
TO = sys.argv[1].split(",")
SUBJECT = sys.argv[2]
TEXT = sys.argv[3]

message = """\
From: %s
To: %s
Subject: %s

%s
""" % (
    FROM,
    ", ".join(TO),
    SUBJECT,
    TEXT,
)

server = smtplib.SMTP("localhost")
server.ehlo()
server.ehlo()  # extra characters to permit edit
server.sendmail(FROM, TO, message)
server.quit()
