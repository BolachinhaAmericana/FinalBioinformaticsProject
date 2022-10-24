FROM python:3.8.2

ADD getSecBdTerm.py .

RUN pip install request biopython

CMD [ "python3","/getSecBdTerm.py" ]