=begin lip

=cut

package Annulus;

use strict;
use warnings; 
use boolean;
use Lip::Pod;
use Math::Trig;
use List::Util qw(max sum);

use Circle;
use constant ORIGIN => [0,0];
use constant PI => 4 * atan2(1, 1);

sub new{
    my ($pkg,$attr)=@_;
    my $self={};
    $self->{link_lengths}=$attr->{link_lengths};
    $self->{link_end_points}=$attr->{link_end_points};
    $self->{center}=ORIGIN;
    ($self->{inner_radius},$self->{outer_radius})=calculate_radii($self->{link_lengths});
    my $inner_circle=new Circle({center=>$self->{center},radius=>$self->{inner_radius}});
    my $outer_circle=new Circle({center=>$self->{center},radius=>$self->{outer_radius}});
    $self->{inner_circle}=$inner_circle;
    $self->{outer_circle}=$outer_circle;
    bless($self,$pkg);
    return $self;
}

=pod 

=head3 calculate_radii()

Using Theorem 8.6.3 we determine the inner and outer radii of the annlus.
The theorem states:
    The reachability region for an n-link arm is an origin-centered annulus
with outer radius r_o=∑_{i=1}^n and inner radius r_i=0 if the longest link 
length l_M is less than or equal to half the total length of the links, and
r_i=l_M-∑_i_{i≠M} l_i otherwise.

=cut

sub calculate_radii{
    my ($link_lengths)=@_;
    my $r_i=0;
    my $r_o=0;
    my $link_sum=sum(@{$link_lengths});
    $r_o=$link_sum;
    my $max_link_length=max(@{$link_lengths});
    unless ($max_link_length <= ($link_sum/2)){
        my @links_no_max=grep {$_ < $max_link_length} @{$link_lengths};
        $r_i=$max_link_length - sum(@links_no_max);
    }
    return ($r_i,$r_o);   
}

sub center{
    my($self)=@_;  
    return $self->{center};
}

sub inner_radius{
    my($self)=@_;  
    return $self->{inner_radius};
}

sub outer_radius{
    my($self)=@_;  
    return $self->{outer_radius};
}

sub link_end_points{
    my($self)=@_; 
    return $self->{link_end_points};
}

=pod

=head3 link_positions()

Returns the signed angle between consecutive links.
The first angle is taken relative to the positive x axis.

=cut

sub link_positions{
    my($self)=@_; 
    my @positions;
    my @link_end_points=@{$self->{link_end_points}};
    my $point_a=shift @link_end_points;#origin
    my $point_b=shift @link_end_points;
    my $angle=atan2($point_b->[1],$point_b->[0]);
    if ($angle < 0){
        $angle=$angle + (2*PI);
    }
    push @positions, $angle;
    my $vertex=$point_b;
    $point_b=shift @link_end_points; 
    while($point_b){
        my $magnitude_a=sqrt(($point_a->[0]-$vertex->[0])**2+($point_a->[1]-$vertex->[1])**2);
        my $magnitude_b=sqrt(($point_b->[0]-$vertex->[0])**2+($point_b->[1]-$vertex->[1])**2);
        my $dot_product=($point_a->[0]-$vertex->[0])*($point_b->[0]-$vertex->[0])+($point_a->[1]-$vertex->[1])*($point_b->[1]-$vertex->[1]);
        $angle=acos($dot_product/($magnitude_a*$magnitude_b));
        if ($angle <= -1*PI){
	        $angle=$angle + (2*PI);
        }
	    elsif ($angle > PI){
	        $angle=$angle - (2*PI);
        }
        push @positions, $angle;
        $point_a=$vertex;
        $vertex=$point_b; 
        $point_b=shift @link_end_points; 
    }
    return \@positions;
}

sub inner_circle{
    my($self)=@_;  
    return $self->{inner_circle};
}

sub outer_circle{
    my($self)=@_;  
    return $self->{outer_circle};
}

sub link_lengths{
    my($self)=@_;
    return $self->{link_lengths};
}

true;

=end lip

=cut