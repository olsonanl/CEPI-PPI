#
# App wrapper for Protein Protein Interface.
# Initial version that does not internally fork and report output_dir; instead
# is designed to be executed by p3x-app-shepherd.
#

use Bio::KBase::AppService::AppScript;
use Bio::KBase::AppService::AppConfig qw(application_backend_dir);

use Clone;
use Cwd qw(abs_path getcwd);
use Data::Dumper;
use feature "switch";
use File::Basename;
use File::Slurp;
use gjoseqlib qw();
use experimental qw(switch);
use IPC::Run qw(run);
use JSON::XS;
use JSON;
use LWP::UserAgent;
use P3DataAPI;
use Path::Tiny qw(path);
use SeedUtils;
use strict;
use URI::Escape;


my $app = Bio::KBase::AppService::AppScript->new(\&run_app, \&preflight);
my $data_api = Bio::KBase::AppService::AppConfig->data_api_url;

$app->run(\@ARGV);

sub run_app
{
    my($app, $app_def, $raw_params, $params) = @_;
    my $begin_time = time();
    print STDERR "Processed parameters for application " . $app->app_definition->{id} . ": ", Dumper($params);
    # print "App-PPI: ", Dumper($app_def, $raw_params, $params);

    ### START Prep ### 
    my $api = P3DataAPI->new();
    my %config_vars;
    #
    # Set up work and staging directories.
    # 
    my $cwd = getcwd();

    my $work_dir = "$cwd/work";
    my $staging_dir = "$cwd/staging";
    my $output_dir = "$cwd/output";
    -d $work_dir or mkdir $work_dir or die "Cannot mkdir $work_dir: $!";
    -d $staging_dir or mkdir $staging_dir or die "Cannot mkdir $staging_dir: $!";
    -d $output_dir or mkdir $output_dir or die "Cannot mkdir $output_dir: $!";
    
    #
    # Write files from user to the staging directory.
    #
    
    # query feature group
    if ($params->{query_feature_group} ne '')
    {
        print($params->{query_feature_group});
        my $ofile = "$staging_dir/query.fasta";
        open(F, ">$ofile") or die "Could not open $ofile";
        my $feature_id_exclude = "";
        # if ($params->{ref_type} eq "feature_id") {
        #     # $feature_id_exclude = $params->{ref_string};
        #     $feature_id_exclude =~ s/^\s+|\s+$//g;
        # }
        # }
        for my $feature_name (@{$params->{query_feature_group}}) {
        my $ids = $api->retrieve_patricids_from_feature_group($feature_name);
        my @ids_new = ();
        for my $id (@$ids){
                # if ($id ne $feature_id_exclude) {
                #     push(@ids_new, $id);
                # }
                push(@ids_new, $id);
            }
            my $seq = "";
                # $seq = $data_api_module->retrieve_protein_feature_sequence(\@ids_new);
                $seq = $api->retrieve_protein_feature_sequence(\@ids_new);
            for my $id (@ids_new) {
                my $out = ">$id\n" . $seq->{$id} . "\n";
                    print F $out;
            $config_vars{query} = $ofile;
            }
        }
    }
    # query fasta file
    if ($params->{query_fasta_file} ne '')
    {
        # load the input fasta file from workspace
        my $out;
        my $local_fasta = "$staging_dir/query.fasta";
        my $ws = Bio::P3::Workspace::WorkspaceClientExt->new;
        if (open(my $out, ">", $local_fasta))
        {

            $ws->copy_files_to_handles(1, undef, [[$params->{query_fasta_file}, $out]]);

            close($out);
        }
        else
        {
            die "Cannot open $local_fasta for writing: $!";
        }
        $config_vars{query} = $out;
    }
    # query keyboard input 
    if ($params->{query_fasta_keyboard_input} ne '')
    {
        my $text_input_file = "$staging_dir/query.fasta";
            open(FH, '>', $text_input_file) or die "Cannot open $text_input_file: $!";
            print FH $params->{query_fasta_keyboard_input};
            close(FH);
            $config_vars{query} = $text_input_file;
    }
    # target feature group
    if ($params->{target_feature_group} ne '')
    {
        my $ofile = "$staging_dir/target.fasta";
        open(F, ">$ofile") or die "Could not open $ofile";
        my $feature_id_exclude = "";
        # if ($params->{ref_type} eq "feature_id") {
        #     # $feature_id_exclude = $params->{ref_string};
        #     $feature_id_exclude =~ s/^\s+|\s+$//g;
        # }
        # }
        for my $feature_name (@{$params->{target_feature_group}}) {
	    # my $ids = $data_api_module->retrieve_patricids_from_feature_group($feature_name);
        my $ids = $api->retrieve_patricids_from_feature_group($feature_name);
        my @ids_new = ();
        for my $id (@$ids){
	        # if ($id ne $feature_id_exclude) {
            #     push(@ids_new, $id);
            # }
            push(@ids_new, $id);
            }
            my $seq = "";
                # $seq = $data_api_module->retrieve_protein_feature_sequence(\@ids_new);
                $seq = $api->retrieve_protein_feature_sequence(\@ids_new);
            for my $id (@ids_new) {
                my $out = ">$id\n" . $seq->{$id} . "\n";
                    print F $out;
            $config_vars{target} = $ofile;
            }
        }
    }
    # target fasta file
    if ($params->{target_fasta_file} ne '')
    {
        # load the input fasta file from workspace
        my $out;
        my $local_fasta = "$staging_dir/target.fasta";
        my $ws = Bio::P3::Workspace::WorkspaceClientExt->new;
        if (open(my $out, ">", $local_fasta))
        {

            $ws->copy_files_to_handles(1, undef, [[$params->{target_fasta_file}, $out]]);

            close($out);
        }
        else
        {
            die "Cannot open $local_fasta for writing: $!";
        }
        $config_vars{target} = $out;
    }
    # target keyboard input
    if ($params->{target_fasta_keyboard_input} ne '')
    {
        my $text_input_file = "$staging_dir/target.fasta";
            open(FH, '>', $text_input_file) or die "Cannot open $text_input_file: $!";
            print FH $params->{target_fasta_keyboard_input};
            close(FH);
            $config_vars{target} = $text_input_file;
    }

    # write remaining values to config in the current working directory

    # $config_vars{cores} = $ENV{P3_ALLOCATED_CPU} // 2;
    # $config_vars{model_path} = $ENV{};
    # hard coded for dev
    $config_vars{params} = $params; 
    # join model_path and pt_model, reformat depending on how you set up the backend directory
    # $config_vars{model_path_future} = application_backend_dir . "CEPI-PPI/models/" . $config_vars{pt_model};
    #$config_vars{model_path} = application_backend_dir . "/CEPI-PPI/models/" . $config_vars{pt_model};
    $config_vars{model_path} = application_backend_dir . "/CEPI-PPI/models/" . "checkpoint-3798";
    $config_vars{threshold} = $params->{threshold};
    $config_vars{input_data_dir} = $staging_dir;
    $config_vars{output_data_dir} = $output_dir;
    $config_vars{work_data_dir} = $work_dir;
    $config_vars{query} = "$staging_dir/query.fasta";
    $config_vars{target} = "$staging_dir/target.fasta";

    my $top = getcwd;
    write_file("$top/config.json", JSON::XS->new->pretty->canonical->encode(\%config_vars));


    # #my @cmd = ("python3", "/nfs/ml_lab/projects/ml_lab/cmann/00_BVBRC_service_development/dev_container/modules/ppi/service-scripts/predict_protein_protein_interface.py");

    # my @cmd = ("python3", "test_script.py");
    # my @cmd = ("python3", "/nfs/ml_lab/projects/ml_lab/cmann/00_BVBRC_service_development/dev_container/modules/ppi/service-scripts/predict_ppi.py", "config.json");
    # my @cmd = ("python3", "/nfs/ml_lab/projects/ml_lab/cmann/00_BVBRC_service_development/dev_container/modules/ppi/service-scripts/test_script.py");
    # my @cmd = ("python3", "/home/nbowers/bvbrc-dev/dev_container/dev/carla_ppi/00_BVBRC_service_development/dev_container/modules/ppi/service-scripts/predict_ppi.py", "config.json");
    
    # production command
    my @cmd = ("predict_ppi", "config.json");
    print STDERR "Run: @cmd\n";
    my $ok = IPC::Run::run(\@cmd);
    if (!$ok)
    {
     die "wrapper command failed $?: @cmd";
    }
    # NB dev  not saving output_dir files
    #save_output_files($app, $top);
}

sub preflight
{
    my($app, $app_def, $raw_params, $params) = @_;

    # ### start From MSA ###
    #     my $query_featuregroups = defined($params->{query_feature_group}) ? $params->{query_feature_group} : undef;
    # if (defined $query_featuregroups) 
    # {
    #     for my $fg (@$query_featuregroups)
    #     {
    #         print "$fg\n";
    #         my $features;
    #         {
    #             $features = $api->retrieve_protein_sequences_from_feature_group($fg);
    #         }
    #         if (defined $features)
    #         {
    #             my $nf = @$features;
    #             $numFeatures = $numFeatures + $nf;
    #         }
    #     }
    # }
    # print("line 65");
    # # print Dumper($query_feature_groups);
    # # # query fasta file
    # my $query_fasta_file = defined($params->{query_fasta_file}) ? $params->{query_fasta_file} : undef;
    # my $totFileSize = 0;
    # if (defined $files) 
    # {
    #     for my $ff (@$files) 
    #     {
    #         print "$ff\n";
    #         my $file_data = $ws->stat($ff->{file});
    #         $totFileSize = $totFileSize + $file_data->size;
    #     }
    # } 
    # lookat the app script  line 36
    # # get number of features from feature groups
    # my $numFeatures = 0;
    # my $featuregroups = defined($params->{feature_groups}) ? $params->{feature_groups} : undef;
    # if (defined $featuregroups) 
    # {
    #     for my $fg (@$featuregroups)
    #     {
    #         print "$fg\n";
    #         my $features;
    #         if ($params->{alphabet} eq "protein")
    #         {
    #             $features = $api->retrieve_protein_sequences_from_feature_group($fg);
    #         } 
    #         else 
    #         {
    #             $features = $api->retrieve_nucleotide_sequences_from_feature_group($fg);
    #         }
    #         if (defined $features)
    #         {
    #             my $nf = @$features;
    #             $numFeatures = $numFeatures + $nf;
    #         }
    #     }
    # }
    # if (exists($params->{feature_list})) {
    #     $numFeatures = $numFeatures + scalar(@{$params->{feature_list}})
    # }

    # # get file sizes
    # my $files = defined($params->{fasta_files}) ? $params->{fasta_files} : undef;
    # my $totFileSize = 0;
    # if (defined $files) 
    # {
    #     for my $ff (@$files) 
    #     {
    #         print "$ff\n";
    #         my $file_data = $ws->stat($ff->{file});
    #         $totFileSize = $totFileSize + $file_data->size;
    #     }
    # } 
        # my $input_type = defined($params->{input_type}) ? $params->{input_type} : undef;

    # my $runtime = 0;
    # my $mem = '';
    # my $mem_threshold = 50000000000; #50GB
    # if (defined $input_type and $input_type eq "input_sequence")
    # {
    #     $runtime = "1800";
    #     $mem = "8GB";
    # }
    # else 
    # {
    #     my $numGroups = $numFeatures + $numGenomes;
    #     if ($numGroups == 0) 
    #     {
    #         $runtime = 3 * 3600;    
    #         $mem = '32GB';
    #     } elsif ($numGroups < 10 and $totFileSize < ($mem_threshold/10)) {
    #         $runtime = 1800; 
    #         $mem = '8GB';
    #     } elsif ($numGroups < 100 and $totFileSize < ($mem_threshold/5))  {
    #         $runtime = 3 * 3600;
    #         $mem = '16GB';
    #     } elsif ($numGroups < 500 and $totFileSize < ($mem_threshold/2)) {
    #         $runtime = 6 * 3600;
    #         $mem = '32GB';
    #     } elsif ($numGroups < 1000 and $totFileSize < ($mem_threshold/2)) {
    #         $runtime = 43200; # 12 hours
    #         $mem = '32GB';
    #     } elsif ($numGroups < 3000 and $totFileSize >= $mem_threshold) {
    #         $runtime = 43200 * 2;
    #         $mem = '64GB';
    #     } else { # <= 5000 genomes
    #         $runtime = 43200 * 4;
    #         $mem = '64GB';
    #     } 
    # }
    # # zero genome ids: default
    # # have no reference for this so just guessing
    # my $pf = {
    #     cpu => 1,
    #     memory => $mem,
    #     runtime => $runtime,
    #     storage => 0,
    #     is_control_task => 0
    # };
    # return $pf;
    ### end FROM MSA ###

}

sub save_output_files
{
    my($app, $output_dir) = @_;
    my %suffix_map = (
        align => 'txt',
    bai => 'bai',
        bam => 'bam',
        csv => 'csv',
        depth => 'txt',
        err => 'txt',
        fasta => "contigs",
        html => 'html',
        out => 'txt',
    png => 'png',
    svg => 'svg',
    tbl => 'tsv',
        tsv => 'tsv',
        txt => 'txt',);
 
    my @suffix_map = map { ("--map-suffix", "$_=$suffix_map{$_}") } keys %suffix_map;
 
    if (opendir(D, $output_dir))
    {
    while (my $p = readdir(D))
    {
        next if ($p =~ /^\./);
        my @cmd = ("p3-cp", "--recursive", @suffix_map, "$output_dir/$p", "ws:" . $app->result_folder);
        print STDERR "saving files to workspace... @cmd\n";
        my $ok = IPC::Run::run(\@cmd);
        if (!$ok)
        {
        warn "Error $? copying output_dir with @cmd\n";
        }
    }
    closedir(D);
    }
}
