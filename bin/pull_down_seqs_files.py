import requests
import pandas as pd
from bs4 import BeautifulSoup
import os
import sys
import pickle
import subprocess

class SeqFilePullDown:
    def __init__(self, dry_run=False):
        print(os.path.abspath(sys.argv[0]))
        self.base_url = "https://www.genoscope.cns.fr/sadc/tarapacific/METABARCODING/18S_V9/18S_V9_1389F_1510R/"
        self.end_url = "https://www.genoscope.cns.fr/sadc/tarapacific/METABARCODING/18S_V9/"
        self.exe_path = os.path.dirname(os.path.abspath(sys.argv[0]))
        self.project_base_dir = os.path.dirname(self.exe_path)
        self.sample_provenance_path = os.path.join(self.project_base_dir, "tara_samples_provenance.csv")
        self.sample_provenance_df = self._make_sample_provenance_df()
        self.list_of_completed_urls = self._get_completed_urls()
        self.authorisation_tup = self._make_auth_tup()
        self.current_url = self.base_url
        self.seq_output_dir = os.path.join(self.project_base_dir, 'seq_files')
        self.dry_run=dry_run
        self.downloaded = self._get_downloaded()
        self.headers = {'User-Agent': 'Benjamin Hume', 'From': 'benjamin.hume@kaust.edu.sa'}

    def _get_completed_urls(self):
        pickle_path = os.path.join(self.exe_path, 'list_of_completed_urls.p')
        if os.path.exists(pickle_path):
            with open(pickle_path, 'rb') as f:
                return pickle.load(f)
        else:
            return []

    def _get_downloaded(self):
        pickle_path = os.path.join(self.exe_path, 'downloaded.p')
        if os.path.exists(pickle_path):
            with open(pickle_path, 'rb') as f:
                return pickle.load(f)
        else:
            return []

    def _pickle_out_object(self, obj, file_name):
        with open(os.path.join(self.exe_path, file_name), 'wb') as f:
            pickle.dump(obj, f)

    def _make_auth_tup(self):
        auth_path = os.path.join(self.exe_path, 'auth.txt')
        with open(auth_path, 'r') as f:
            auth_lines = [line.rstrip() for line in f]
        return (auth_lines[0], auth_lines[1])

    def _make_sample_provenance_df(self):
        # The SAMPLE ID is in the sequencing file name so we will be able to use this to get the latitute
        df = pd.read_csv(self.sample_provenance_path)
        df.set_index(keys='SAMPLE ID', drop=True, inplace=True)
        df.rename(columns={'EVENT latitude start (North)':'lat', 'EVENT longitude start (East)': 'lon'}, inplace=True)
        return df

    def download_coral_files(self):
        """Recursively go through the file structure down from the base url
        When you find a fastq.gz, extract the sample name, check it in the sample_provenance_df
        and if it is coral, down load it. Once we have been through all viable links, we go one directory up and
        continue the search. If we are in the base dir, and there are no further links to visit then
        it will be time to call it quits."""
        self._walk()

    def _walk(self):
        while True:
            print(f'Current directory: {self.current_url}')
            soup = BeautifulSoup(requests.get(self.current_url, auth=self.authorisation_tup, headers=self.headers).text, features="html.parser")
            # Strip out the table links and those that have already been visited
            links_to_visit = [link.string for link in soup.find_all('a') if
                              ((link.string not in ['Name', 'Last modified', 'Size', 'Description', 'Parent Directory', 'NEGATIVE_CONTROLS/']) and
                               ('/'.join([self.current_url.strip('/'), link.string]) not in self.list_of_completed_urls))]
            fastq_gz_list = [link.string for link in links_to_visit if 'fastq.gz' in link.string]
            links_to_visit = [link for link in links_to_visit if link not in fastq_gz_list]
            if fastq_gz_list:
                self._assess_fastq_gz_files(fastq_gz_list)

            # There are no fastq_gz files so we need to coninue walking if there are any valid links
            if links_to_visit:
                self.current_url = '/'.join([self.current_url.strip('/'), links_to_visit[0]])
                continue
            else:
                # If there are no further links to visit for the current url
                # then we add it to the list_of_completed_urls
                self.list_of_completed_urls.append(self.current_url)
                self._pickle_out_object(self.list_of_completed_urls, 'list_of_completed_urls.p')

            # If we get here then we didn't find a link to visit and it is time to head up a directory
            # We will need to check to see that we aren't back at the base directory
            # Go up one directory unless we are already back at the root directory
            print('All links visited, walking up')
            self.current_url = '/'.join(self.current_url.split('/')[:-2]) + '/'
            if self.current_url in self.end_url:
                # We have completed the walk and it is time to exit out
                print(f'{len(self.list_of_completed_urls)} URLs visited')
                print(f'{len(self.downloaded)} files downloaded')
                break
            else:
                # restart the walk
                continue

    def _assess_fastq_gz_files(self, fastq_gz_list):
        for fastq_file in fastq_gz_list:
            print(f'Checking if {fastq_file} is coral')
            # Check to see that it is not in the already visited list
            # Also check to see that it has not already been downloaded
            # We will have written out a simple text file to shown that the download completed. Look for this

            download_url = '/'.join([self.current_url.strip('/'), fastq_file])
            if download_url not in self.list_of_completed_urls:
                if not os.path.exists(os.path.join(self.seq_output_dir, fastq_file.replace('.fastq.gz', '.txt'))):
                    # Then we need to check it to see if it is type coral
                    if self._if_sample_is_coral(fastq_file):
                        if self.dry_run:
                            print(f'Downloading {fastq_file}')
                            self.downloaded.append(fastq_file)
                            self.list_of_completed_urls.append(download_url)
                            print(f'Download complete')
                        else:
                            # Then we want to download this file
                            # For some reason the files are downloading uncompressed
                            # We will compress them at a later date as we can multiprocess this.
                            print(f'Downloading {fastq_file}')
                            r = requests.get(download_url, auth=self.authorisation_tup)
                            fastq_path = os.path.join(self.seq_output_dir, fastq_file.replace('.gz', ''))
                            open(fastq_path, 'wb').write(r.content)
                            # subprocess.run(['gzip', fastq_path])
                            open(
                                os.path.join(self.seq_output_dir, fastq_file.replace('.fastq.gz', '.txt')), 'w'
                            ).write('complete\n')
                            # Then add the url to the list_of_visited_urls
                            self.list_of_completed_urls.append(download_url)
                            self._pickle_out_object(self.list_of_completed_urls, 'list_of_completed_urls.p')
                            self.downloaded.append(fastq_file)
                            self._pickle_out_object(self.downloaded, 'downloaded.p')
                            print(f'Download complete')
                    else:
                        # This is not a coral sample and we don't need to download it
                        print(f'{fastq_file} is not a coral sample')
                        self.list_of_completed_urls.append(download_url)
                        self._pickle_out_object(self.list_of_completed_urls, 'list_of_completed_urls.p')
                else:
                    print(f'{fastq_file} already exists locally')
                    self.list_of_completed_urls.append(download_url)
                    self._pickle_out_object(self.list_of_completed_urls, 'list_of_completed_urls.p')

    def _if_sample_is_coral(self, fastq_name):
        sample_name = '_'.join(fastq_name.split('_')[:2])
        sample_type = self.sample_provenance_df.at[sample_name, 'SAMPLE MATERIAL LABEL']
        print(f'{fastq_name} is {sample_type}')
        return sample_type == 'CORAL'

sfpd = SeqFilePullDown(dry_run=False)
sfpd.download_coral_files()