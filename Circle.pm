=begin lip

=cut

package Circle;

use strict;
use warnings; 
use boolean;
use Lip::Pod;

sub new{
    my ($pkg,$attr)=@_;
    my $self={};
    $self->{center}=$attr->{center};
    $self->{radius}=$attr->{radius};
    bless($self,$pkg);
    return $self;
}

sub center{
    my($self)=@_;  
    return $self->{center};
}

sub radius{
    my($self)=@_;  
    return $self->{radius};
}

true;

=end lip

=cut