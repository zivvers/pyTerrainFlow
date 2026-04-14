import tomllib
from dataclasses import dataclass
from pathlib import Path

@dataclass(frozen=True)
class Config:
    input_dir: Path
    base_dir: Path
    output_dir: Path
    #road_gpkg: Path
    #major_road_gpkg: Path
    #oregon_border: Path


def load_config(config_path = "config.toml") -> Config:
    config_path = Path(config_path)

    with config_path.open("rb") as f:
        raw = tomllib.load(f)
    
    base_dir = Path( raw["paths"]["base_dir"] ).resolve()

    input_dir = base_dir / raw["paths"]["input_dir"]

    output_dir = base_dir / raw["paths"]["output_dir"]

    return Config( base_dir=base_dir
        , output_dir=output_dir
        , input_dir=input_dir )
