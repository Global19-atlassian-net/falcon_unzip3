import logging
import traceback
import subprocess
from time import gmtime, strftime


def log(message, logger):
  """Return the current date and time.
  """
  timestamp = strftime("%Y/%m/%d %H:%M:%S", gmtime());
  msg = '[{}] {}'.format(timestamp, message)
  logger.log(logging.WARNING-1, msg)


def execute_command(command, logger, dry_run=False):
  """
  Executes a given shell command and logs the attempt.
  """
  if dry_run:
      log('Executing "{}" (dryrun={})'.format(command, dry_run), logger)
      logger.info('')
      return 0
  log('Executing "{}"'.format(command), logger)
  rc = subprocess.call(command, shell=True, encoding='ascii')
  if (rc != 0):
      msg = ' subprocess call returned error code: {}.'.format(rc)
      log(msg, logger)
      logger.exception('{}<-"{}"'.format(rc, command))
      raise Exception(msg)
  log(' Finished subprocess.', logger)
  return rc


def execute_command_with_ret(command, logger, dry_run=False):
  """
  Executes a given shell command and logs the attempt, but also gets
  the stdout, stderr and the return code of the program.
  """
  log('Executing "{}" (dryrun={}).'.format(command, dry_run), logger)
  if dry_run:
    return [0, '', '']
  p = subprocess.Popen(command, shell=True, encoding='ascii', stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  [output, err] = p.communicate()
  rc = p.returncode

  log(' Finished subprocess, code={}'.format(rc), logger)

  return [rc, output, err]
