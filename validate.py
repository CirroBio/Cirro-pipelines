import json
from json.decoder import JSONDecodeError
import logging
import jsonschema
from pathlib import Path

logging.basicConfig(format='[validate] - %(levelname)s - %(message)s',
                    level=logging.INFO)

logger = logging.getLogger()

schemas = {}


def load_schemas():
    for file in list(Path('schemas').glob('*.schema.json')):
        with file.open() as f:
            schema_name = file.name.replace('.schema.json', '')
            schemas[schema_name] = json.load(f)


def load_items(file_name):
    items = {}
    files = list(Path('.').glob(f'*/**/{file_name}'))
    for file in files:
        with file.open() as f:
            try:
                items[file] = json.load(f)
            except JSONDecodeError as err:
                logger.error(f"Error loading {file}:")
                raise err
    return items


def validate_items(schema_name):
    file_name = f'{schema_name}.json'
    items_to_validate = load_items(file_name)
    has_errors = False
    for name, contents in items_to_validate.items():
        try:
            jsonschema.validate(
                instance=contents,
                schema=schemas.get(schema_name)
            )
        except jsonschema.exceptions.ValidationError as err:
            has_errors = True
            logger.warning(
                f"{name} JSON doesn't validate against '{schema_name}' schema, {err}" # noqa
            )
    if has_errors:
        raise RuntimeError("JSON Schema errors found, see log")
    else:
        logger.info(f"All {schema_name} objects pass JSON validation")


def validate_process_records():
    process_records = load_items('process-definition.json')
    process_ids = [process['id'] for process in process_records.values()]
    has_errors = False

    # Check if process IDs are unique
    duplicates = set([x for x in process_ids if process_ids.count(x) > 1])
    if len(duplicates):
        has_errors = True
        logger.warning(f"Duplicate process IDs found {duplicates}")

    for process in process_records.values():
        # Check if child process IDs are valid
        for child_process in process['childProcessIds']:
            if child_process not in process_ids:
                has_errors = True
                logger.warning(f"child process {child_process} in {process['id']} is not valid")

        # Check file extensions
        for field in ["paramMapJson", "formJson", "webOptimizationJson"]:
            if process[field] and not process[field].endswith('.json'):
                has_errors = True
                logger.warning(f"Invalid path {process[field]} in {process['id']}")

    if has_errors:
        raise RuntimeError("Process record errors found, see log")
    else:
        logger.info("All processes pass validation")


def run_validation():
    load_schemas()

    validate_items('process-form')
    validate_items('process-definition')
    validate_items('process-input')
    validate_process_records()


if __name__ == '__main__':
    run_validation()
