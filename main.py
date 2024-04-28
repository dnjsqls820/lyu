# from imp import reload
# from pyparsing import Char
from token import OP
from fastapi import FastAPI, Depends, Request
import pandas as pd
import sqlalchemy as sa
from dataclasses import dataclass, Field
from sqlalchemy.orm import sessionmaker
from sqlalchemy import VARCHAR, Column, Integer, String, func, or_, and_
import json
from sqlalchemy.orm import Session

from typing import List, Text, Optional
from pydantic import BaseModel
from sqlalchemy.orm import declarative_base

# from sqlalchemy.ext.declarative import declarative_base
import uvicorn
from fastapi.responses import JSONResponse
from fastapi_pagination import Page, Params, paginate

Base = declarative_base()
# Local_DB = "mysql+pymysql://root@127.0.0.1:3306/new_schema"
Local_DB = "mysql+pymysql://root:isa4986@166.104.112.52:3306/immunoresource?charset=utf8"
engine = sa.create_engine(Local_DB)
app = FastAPI()
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)


def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()


class t_MAP_statistics(Base):
    __tablename__ = "MAP_statistics"
    id = Column(Integer, primary_key=True, index=True)
    tissue_name = Column(String(64), unique=True)
    sample_count = Column(Integer)
    MAP_count_c = Column(Integer)
    MAP_count_nc = Column(Integer)


class r_MAP_statistics(BaseModel):
    id: int
    tissue_name: str
    sample_count: int
    MAP_count_c: int
    MAP_count_nc: int


class t_total_statistics(Base):
    __tablename__ = "total_statistics"
    total_samples = Column(Integer, primary_key=True)
    total_tissues = Column(Integer)
    total_cancer_types = Column(Integer)
    total_peptides_c = Column(Integer)
    total_peptides_nc = Column(Integer)


class r_total_statistics(BaseModel):
    total_samples: int
    total_tissues: int
    total_cancer_types: int
    total_peptides_c: int
    total_peptides_nc: int


class Psm(Base):
    __tablename__ = "PSM_table"
    id = Column(Integer, primary_key=True)
    gene_name = Column(VARCHAR(1000))
    peptide_sequence = Column(VARCHAR(50))
    hla_type = Column(VARCHAR(50))
    event = Column(VARCHAR(50))
    tissue = Column(VARCHAR(50))
    disease = Column(VARCHAR(50))
    peptide_type = Column(VARCHAR(20))
    hla_allele = Column(VARCHAR(200))


@app.get("/statistic-info")
def get_data(db: Session = Depends(get_db)):
    total_data = db.query(t_total_statistics).first()
    MAP_data = db.query(t_MAP_statistics).all()
    result = {"total_statistics": total_data, "MAP_statistics": MAP_data}
    return result


class Item(BaseModel):
    page: int
    per_page: int
    tissue_exclusive: str
    allele_exclusive: str
    peptide_type: str
    hla_type: str
    column_name: Optional[str] = None
    sort: Optional[str] = None
    gene_list: Optional[list] = None
    sequence_list: Optional[list] = None
    tissue_type_list: Optional[list] = None
    event_type_list: Optional[list] = None
    hla_allele_list: Optional[list] = None


@app.post("/peptide-table/")
def test(request: Item, db: Session = Depends(get_db)):

    offset = (request.page - 1) * request.per_page
    query = db.query(Psm)
    if request.tissue_exclusive == "True":
        query = query.filter(func.length(Psm.tissue) - func.length(func.replace(Psm.tissue, "|", "")) + 1 == 1)
    if request.allele_exclusive == "True":
        query = query.filter(func.length(Psm.hla_allele) - func.length(func.replace(Psm.hla_allele, "|", "")) + 1 == 1)
    if request.peptide_type != "all":
        query = query.filter(Psm.peptide_type == request.peptide_type)
    if request.hla_type != "all":
        query = query.filter(Psm.hla_type == f"HLA-Class {request.hla_type}")

    if request.gene_list:
        conditions = [Psm.gene_name.like(f"%{value}%") for value in request.gene_list]
        query = query.filter(or_(*conditions))

    if request.sequence_list:
        conditions = [Psm.peptide_sequence.like(f"%{value}%") for value in request.sequence_list]
        query = query.filter(or_(*conditions))

    if request.tissue_type_list:
        conditions = [Psm.tissue.like(f"%{value}%") for value in request.tissue_type_list]
        query = query.filter(or_(*conditions))
    if request.event_type_list:
        conditions = [Psm.event.like(f"%{value}%") for value in request.event_type_list]
        query = query.filter(or_(*conditions))
    if request.hla_allele_list:
        conditions = [Psm.hla_allele.like(f"%{value}%") for value in request.hla_allele_list]
        query = query.filter(or_(*conditions))

    if request.sort and request.column_name:
        column = getattr(Psm, request.column_name)
        if request.sort == "desc":
            query = query.order_by(column.desc())
        if request.sort == "asc":
            query = query.order_by(column.asc())

    total_count = db.scalar(db.query(func.count(Psm.id)))

    # 전체 페이지
    query = query.offset(offset).limit(request.per_page)
    data = query.all()
    search_cnt = len(data)
    total_pages = round(search_cnt / request.per_page)
    pages = 1 if 1 > total_pages else total_pages
    pagination = {
        "peptide_table_data": data,
        "pagination": {
            "page": request.page,
            "pages": pages,
            "prev_page": "null" if request.page == 1 else request.page - 1,
            "next_page": "null" if request.page == pages else request.page + 1,
            "total": total_count,
        },
    }

    return pagination


if __name__ == "__main__":

    uvicorn.run("main:app", host="0.0.0.0", port=8000, reload=True)
