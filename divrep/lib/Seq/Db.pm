package Seq::Db;
use strict;
use warnings;

use DBI;
BEGIN {
    $ENV{DBIC_DONT_VALIDATE_RELS} = 1;
}

sub new {
    my($class, $type, $recreate) = @_;
    my($dbname, $dbuser, $dbpass)
            = ($type->dbname, $type->dbuser, $type->dbpass);
    my $dsn = sprintf(
        'DBI:%s:database=%s;host=%s',
        'mysql', $dbname // die("No db specified"), 'localhost',
    );
    my $db = Seq::Db::Schema->connect($dsn, $dbuser, $dbpass)
            // die "Schema failed: $DBI::errstr";
    $db->deploy($recreate ? { add_drop_table => 1 } : ());
    return bless {
        type => $type,
        db => $db,
    }, $class;
}

sub db { shift->{db} }
sub resultset { shift->db->resultset(@_) }
sub type { shift->{type} }

package Seq::Db::Schema {
    use Seq::TauG;
    use Seq::TauF;
    use Seq::Run;
    use parent qw{ DBIx::Class::Schema };
    __PACKAGE__->load_classes({
        Seq => [qw{ TauG TauF Run }],
    });
    sub deployment_statements {
        my($self, @args) = @_;
        my @results = $self->next::method(@args);
        s{CREATE TABLE}{CREATE TABLE IF NOT EXISTS}g for @results;
        return @results;
    }
};

1;
