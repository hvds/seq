package Seq::Db;
use strict;
use warnings;

use DBI;
BEGIN {
    $ENV{DBIC_DONT_VALIDATE_RELS} = 1;
}

sub new {
    my($class, $args) = @_;
    my $dsn = sprintf(
        'DBI:%s:database=%s;host=%s',
        'mysql', $args->{db} // die("No db specified"), 'localhost',
    );
    my $db = Seq::Db::Schema->connect($dsn, @$args{qw{ user passwd }})
            // die "Schema failed: $DBI::errstr";
    $db->deploy($args->{recreate} ? { add_drop_table => 1 } : ());
    return bless {
        args => $args,
        db => $db,
    }, $class;
}

sub db { shift->{db} }
sub resultset { shift->db->resultset(@_) }
sub gprio {
    my($self, $count) = @_;
    return Seq::TauG->prio($self, $count);
}

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
